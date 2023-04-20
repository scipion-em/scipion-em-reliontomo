# *
# * Authors:     Scipion Team
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307  USA
# *
# *  All comments concerning this program package may be sent to the
# *  e-mail address 'scipion-users@lists.sourceforge.net'
# *
# **************************************************************************

import logging
logger = logging.getLogger(__name__)
from os import symlink
from emtable import Table
from pwem.emlib.image import ImageHandler
import pyworkflow.utils as pwutils
from pwem.objects import Transform
from pyworkflow.object import Float
from pyworkflow.utils import removeBaseExt, getParentFolder
from reliontomo.constants import FILE_NOT_FOUND, COORD_X, COORD_Y, COORD_Z, SUBTOMO_NAME, \
    TILT_PRIOR, PSI_PRIOR, TOMO_NAME_30, CTF_MISSING_WEDGE, \
    CLASS_NUMBER, MRC, PARTICLES_TABLE
from reliontomo.convert.convertBase import WriterTomo, getTransformInfoFromCoordOrSubtomo, ReaderTomo, \
    getTransformMatrixFromRow
from reliontomo.utils import getAbsPath, _gen2LevelBaseName
from scipion.install.funcs import mkdir
from os.path import join
from tomo.constants import BOTTOM_LEFT_CORNER
from tomo.objects import SubTomogram, Coordinate3D, TomoAcquisition


class Writer(WriterTomo):
    """ Helper class to convert from Scipion SetOfImages subclasses
    with star file format previous to Relion>3.1, but providing the same
     interface as the new Writer class.
    """

    def __init__(self,  **kwargs):
        super().__init__(**kwargs)

    def subtomograms2Star(self, subtomoSet, subtomosStar):

        logger.info("Converting scipion subtomograms to relion 3 format.")
        currentTomo = ''
        ih = ImageHandler()
        tomoTable = Table(columns=self.starHeaders)
        tmpDir = pwutils.getParentFolder(subtomosStar)
        for subtomo in subtomoSet.iterSubtomos():
            if pwutils.getExt(subtomo.getFileName().replace(':' + MRC, '')) != '.' + MRC:
                mrcDir = join(tmpDir, pwutils.removeBaseExt(subtomo.getVolName()))
                if currentTomo != subtomo.getVolName():
                    mkdir(mrcDir)
                mrcFile = join(mrcDir, pwutils.replaceBaseExt(subtomo.getFileName(), MRC))
                ih.convert(subtomo.getFileName(), mrcFile)
            angles, shifts = getTransformInfoFromCoordOrSubtomo(subtomo, 1)
            magn = subtomo.getAcquisition().getMagnification()

            rlnCoordinateX = 0
            rlnCoordinateY = 0
            rlnCoordinateZ = 0
            coord3D = subtomo.getCoordinate3D()

            ctfFile = getattr(coord3D, '_3dcftMrcFile', None)
            if ctfFile:
                ctfFile = ctfFile.get()

            rlnMicrographName = subtomo.getVolName()
            if coord3D:
                rlnCoordinateX = coord3D.getX(BOTTOM_LEFT_CORNER)
                rlnCoordinateY = coord3D.getY(BOTTOM_LEFT_CORNER)
                rlnCoordinateZ = coord3D.getZ(BOTTOM_LEFT_CORNER)

            rlnImageName = subtomo.getFileName().replace(':' + MRC, '')
            rlnCtfImage = ctfFile if ctfFile else FILE_NOT_FOUND
            rlnMagnification = magn if magn else 10000 #64000
            rlnDetectorPixelSize = subtomo.getSamplingRate()
            rlnAngleRot = angles[0]
            rlnAngleTilt = angles[1]
            rlnAnglePsi = angles[2]
            rlnOriginX = shifts[0]
            rlnOriginY = shifts[1]
            rlnOriginZ = shifts[2]
            rlnTiltPrior = subtomo._tiltPriorAngle.get() if hasattr(subtomo, '_tiltPriorAngle') else rlnAngleTilt
            rlnPsiPrior = subtomo._psiPriorAngle.get() if hasattr(subtomo, '_psiPriorAngle') else rlnAnglePsi
            # Add row to the table which will be used to generate the STAR file
            fieldsToAdd = [rlnMicrographName,
                           rlnCoordinateX,
                           rlnCoordinateY,
                           rlnCoordinateZ,
                           rlnImageName,
                           rlnCtfImage,
                           rlnMagnification,
                           rlnDetectorPixelSize,
                           rlnAngleRot,
                           rlnAngleTilt,
                           rlnTiltPrior,
                           rlnAnglePsi,
                           rlnPsiPrior,
                           rlnOriginX,
                           rlnOriginY,
                           rlnOriginZ]

            tomoTable.addRow(*fieldsToAdd)

        # Write the STAR file
        tomoTable.write(subtomosStar)


class Reader(ReaderTomo):

    def __init__(self, starFile, dataTable, **kwargs):
        super().__init__(starFile, dataTable)

    @staticmethod
    def gen3dCoordFromStarRow(row, precedentDict, scaleFactor=1):
        coordinate3d = Coordinate3D()
        tomoNameFromStar = row.get(TOMO_NAME_30)
        x = float(row.get(COORD_X, 0))
        y = float(row.get(COORD_Y, 0))
        z = float(row.get(COORD_Z, 0))
        tomo = precedentDict[removeBaseExt(tomoNameFromStar)]
        coordinate3d.setVolume(tomo)
        ctf3d = row.get(CTF_MISSING_WEDGE, FILE_NOT_FOUND)
        coordinate3d.setX(x * scaleFactor, BOTTOM_LEFT_CORNER)
        coordinate3d.setY(y * scaleFactor, BOTTOM_LEFT_CORNER)
        coordinate3d.setZ(z * scaleFactor, BOTTOM_LEFT_CORNER)
        coordinate3d._3dcftMrcFile = ctf3d  # Used for the ctf3d generation in Relion
        coordinate3d.setMatrix(getTransformMatrixFromRow(row))

        return coordinate3d

    def starFile2Coords3D(self, coordsSet, precedentsSet, scaleFactor):
        precedentDict = {removeBaseExt(tomo.getFileName()): tomo.clone() for tomo in precedentsSet}
        for row in self.dataTable:
            coordsSet.append(self.gen3dCoordFromStarRow(row, precedentDict, scaleFactor))

    def starFile2SubtomogramsImport(self, subtomoSet, coordSet, linkedSubtomosDir, starFilePath):
        samplingRate = subtomoSet.getSamplingRate()
        for row, coordinate3d in zip(self.dataTable, coordSet):
            subtomo = SubTomogram()
            transform = Transform()
            origin = Transform()

            # Files
            tomoName = row.get(TOMO_NAME_30, FILE_NOT_FOUND)
            subtomoName = row.get(SUBTOMO_NAME, FILE_NOT_FOUND)
            linkedSubtomoName = join(linkedSubtomosDir, _gen2LevelBaseName(subtomoName))
            symlink(getAbsPath(starFilePath, subtomoName), linkedSubtomoName)  # Link the subtomos to the extra folder

            # Subtomograms
            tiltPrior = row.get(TILT_PRIOR, 0)
            psiPrior = row.get(PSI_PRIOR, 0)
            transform.setMatrix(coordinate3d.getMatrix())

            subtomo.setVolName(tomoName)
            subtomo.setFileName(linkedSubtomoName)
            subtomo.setCoordinate3D(coordinate3d)
            subtomo.setVolId(coordinate3d.getVolId())
            subtomo.setTransform(transform)
            subtomo.setAcquisition(TomoAcquisition())
            subtomo.setClassId(row.get(CLASS_NUMBER, 0))
            subtomo.setSamplingRate(samplingRate)
            subtomo._tiltPriorAngle = Float(tiltPrior)
            subtomo._psiPriorAngle = Float(psiPrior)
            subtomo.setOrigin(origin)

            # Add current subtomogram to the output set
            subtomoSet.append(subtomo)

        # Set the set of coordinates which corresponds to the current set of subtomograms
        subtomoSet.setCoordinates3D(coordSet)



