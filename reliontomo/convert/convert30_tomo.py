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
from os import symlink
from emtable import Table
from pwem.convert import transformations
from pwem.emlib.image import ImageHandler
import pyworkflow.utils as pwutils
from pwem.objects import Transform
from pyworkflow.object import Float
from pyworkflow.utils import removeBaseExt
from relion.convert.convert_base import WriterBase, ReaderBase
from reliontomo.constants import FILE_NOT_FOUND, COORD_X, COORD_Y, COORD_Z, SUBTOMO_NAME, ROT, TILT, \
    TILT_PRIOR, PSI, PSI_PRIOR, SHIFTX, SHIFTY, SHIFTZ, TOMO_NAME_30, CTF_MISSING_WEDGE, MAGNIFICATION, PIXEL_SIZE, \
    CLASS_NUMBER
from reliontomo.utils import getAbsPath, _gen2LevelBaseName, getTransformMatrix
from scipion.install.funcs import mkdir
import numpy as np
from os.path import join
from pwem.convert.transformations import translation_from_matrix, euler_from_matrix
from tomo.constants import BOTTOM_LEFT_CORNER
from tomo.objects import SubTomogram, Coordinate3D, TomoAcquisition

RELION_TOMO_LABELS = [TOMO_NAME_30,
                      COORD_X,
                      COORD_Y,
                      COORD_Z,
                      SUBTOMO_NAME,
                      CTF_MISSING_WEDGE,
                      MAGNIFICATION,
                      PIXEL_SIZE,
                      ROT,
                      TILT,
                      TILT_PRIOR,
                      PSI,
                      PSI_PRIOR,
                      SHIFTX,
                      SHIFTY,
                      SHIFTZ]


class Writer(WriterBase):
    """ Helper class to convert from Scipion SetOfImages subclasses
    with star file format previous to Relion>3.1, but providing the same
     interface as the new Writer class.
    """

    @staticmethod
    def writeSetOfSubtomograms(subtomoSet, subtomosStar, isPyseg=False):
        currentTomo = ''
        MRC = 'mrc'
        ih = ImageHandler()
        tomoTable = _createStarTomoTable(isPyseg)
        tmpDir = pwutils.getParentFolder(subtomosStar)
        for subtomo in subtomoSet:
            if pwutils.getExt(subtomo.getFileName().replace(':' + MRC, '')) != '.' + MRC:
                mrcDir = join(tmpDir, pwutils.removeBaseExt(subtomo.getVolName()))
                if currentTomo != subtomo.getVolName():
                    mkdir(mrcDir)
                mrcFile = join(mrcDir, pwutils.replaceBaseExt(subtomo.getFileName(), MRC))
                ih.convert(subtomo.getFileName(), mrcFile)
            angles, shifts = _getTransformInfoFromSubtomo(subtomo)
            magn = subtomo.getAcquisition().getMagnification()
            ctfFile = getattr(subtomo.getCoordinate3D(), '_3dcftMrcFile', None)
            if ctfFile:
                ctfFile = ctfFile.get()
            rlnMicrographName = subtomo.getVolName()
            rlnCoordinateX = subtomo.getCoordinate3D().getX(BOTTOM_LEFT_CORNER)
            rlnCoordinateY = subtomo.getCoordinate3D().getY(BOTTOM_LEFT_CORNER)
            rlnCoordinateZ = subtomo.getCoordinate3D().getZ(BOTTOM_LEFT_CORNER)
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
            rlnTiltPsi = subtomo._psiPriorAngle.get() if hasattr(subtomo, '_psiPriorAngle') else rlnAnglePsi
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
                           rlnTiltPsi,
                           rlnOriginX,
                           rlnOriginY,
                           rlnOriginZ]
            if isPyseg:
                fieldsToAdd = [rlnMicrographName,
                               rlnCoordinateX,
                               rlnCoordinateY,
                               rlnCoordinateZ,
                               rlnImageName,
                               rlnCtfImage,
                               rlnAngleRot,
                               rlnAngleTilt,
                               rlnAnglePsi,
                               rlnOriginX,
                               rlnOriginY,
                               rlnOriginZ]

            tomoTable.addRow(*fieldsToAdd)

        # Write the STAR file
        tomoTable.write(subtomosStar)


class Reader(ReaderBase):

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.dataTable = Table()

    def read(self, starFile):
        self.dataTable.read(starFile)

    @staticmethod
    def gen3dCoordFromStarRow(row, precedentsSet, precedentIdList, scaleFactor):
        coordinate3d = Coordinate3D()
        tomoId = removeBaseExt(row.get(TOMO_NAME_30))
        x = float(row.get(COORD_X, 0))
        y = float(row.get(COORD_Y, 0))
        z = float(row.get(COORD_Z, 0))
        volId = precedentIdList.index(tomoId) + 1  # Set indices begin in 1
        coordinate3d.setVolume(precedentsSet[volId])
        coordinate3d.setVolId(volId)
        ctf3d = row.get(CTF_MISSING_WEDGE, FILE_NOT_FOUND)
        coordinate3d.setTomoId(tomoId)
        coordinate3d.setX(x * scaleFactor, BOTTOM_LEFT_CORNER)
        coordinate3d.setY(y * scaleFactor, BOTTOM_LEFT_CORNER)
        coordinate3d.setZ(z * scaleFactor, BOTTOM_LEFT_CORNER)
        coordinate3d._3dcftMrcFile = ctf3d  # Used for the ctf3d generation in Relion
        coordinate3d.setMatrix(getTransformMatrixFromRow(row))

        return coordinate3d

    def starFile2Coords3D(self, coordsSet, precedentsSet, scaleFactor):
        precedentIdList = [tomo.getTsId() for tomo in precedentsSet]
        for row in self.dataTable:
            coordsSet.append(self.gen3dCoordFromStarRow(row, precedentsSet, precedentIdList, scaleFactor))

    def starFile2Subtomograms(self, subtomoSet, coordSet, linkedSubtomosDir, starFilePath):
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


def _createStarTomoTable(isPyseg):

    cols = RELION_TOMO_LABELS
    # Pyseg post-rec only works if the magnification, pixel size and the prior angles aren't
    # present in the star file
    if isPyseg:
        cols = [TOMO_NAME_30,
                COORD_X,
                COORD_Y,
                COORD_Z,
                SUBTOMO_NAME,
                CTF_MISSING_WEDGE,
                ROT,
                TILT,
                PSI,
                SHIFTX,
                SHIFTY,
                SHIFTZ]
    return Table(columns=cols)


def _getCTFFileFromSubtomo(subtomo):
    try:
        return subtomo.getCoordinate3D()._3dcftMrcFile.get()
    except:
        return 'Unavailable'


def _getTransformInfoFromSubtomo(subtomo, calcInv=True):
    angles = [0, 0, 0]
    shifts = [0, 0, 0]
    T = subtomo.getTransform()

    if T:  # Alignment performed before
        M = subtomo.getTransform().getMatrix()
        shifts = translation_from_matrix(M)
        if calcInv:
            shifts = -shifts
            M = np.linalg.inv(M)

        angles = -np.rad2deg(euler_from_matrix(M, axes='szyz'))

    return angles, shifts


def getTransformMatrixFromRow(row, invert=True):
    shiftx = row.get(SHIFTX, 0)
    shifty = row.get(SHIFTY, 0)
    shiftz = row.get(SHIFTZ, 0)
    tilt = row.get(TILT, 0)
    psi = row.get(PSI, 0)
    rot = row.get(ROT, 0)

    return getTransformMatrix(shiftx, shifty, shiftz, rot, tilt, psi, invert)

