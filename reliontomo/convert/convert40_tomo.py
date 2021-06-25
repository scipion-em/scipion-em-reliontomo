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

from pwem.emlib.image import ImageHandler
import pyworkflow.utils as pwutils
from relion.convert.convert_base import WriterBase
from reliontomo.constants import TOMO_NAME, TILT_SERIES_NAME, CTFPLOTTER_FILE, IMOD_DIR, FRACTIONAL_DOSE, \
    ACQ_ORDER_FILE, CULLED_FILE, CTF_IMAGE, SUBTOMO_NAME, COORD_X, COORD_Y, COORD_Z, SHIFTX, SHIFTY, SHIFTZ, ROT, \
    TILT, PSI, TILT_PRIOR, PSI_PRIOR, CLASS_NUMBER
from scipion.install.funcs import mkdir
import numpy as np
from os.path import abspath, join
from pwem.convert.transformations import translation_from_matrix, euler_from_matrix
from relion.convert import Table

# Star file fields
from tomo.constants import BOTTOM_LEFT_CORNER


class Writer(WriterBase):
    """ Helper class to convert from Scipion SetOfTomograms and SetOfSubTomograms to star files ."""

    def writeSetOfTomograms(self, tomoSet, tsSet, ctfPlotterDir, eTomoDir, outStarFileName):
        tomoTable = Table(columns=self._getTomogramStarFileLabels())
        # rr, tt = zip(*[(i * 10, i * 12) for i in xrange(4)])
        tsFileList, tsIdFromTsSet = zip(*[(ts.getFileName(), ts.getTsId()) for ts in tsSet])
        for tomo in tomoSet:
            tsId = tomo.getTsId()
            tomoName = tomo.getFileName()
            tsName = tsFileList.index(tsId)
            ctfPlotterFile = self._getCtfPlotterFile(tsId, ctfPlotterDir)
            eTomoTsSubDir = join(eTomoDir, tsId)
            # fractionalDose =


    def writeSetOfSubtomograms(self, subtomoSet, subtomosStar, isPyseg=False, **kwargs):
        currentTomo = ''
        MRC = 'mrc'
        ih = ImageHandler()
        tomoTable = self._createStarTomoTable(isPyseg)
        tmpDir = pwutils.getParentFolder(subtomosStar)
        for subtomo in subtomoSet:
            if pwutils.getExt(subtomo.getFileName().replace(':' + MRC, '')) != '.' + MRC:
                mrcDir = join(tmpDir, pwutils.removeBaseExt(subtomo.getVolName()))
                if currentTomo != subtomo.getVolName():
                    mkdir(mrcDir)
                mrcFile = join(mrcDir, pwutils.replaceBaseExt(subtomo.getFileName(), MRC))
                ih.convert(subtomo.getFileName(), mrcFile)
            angles, shifts = self._getTransformInfoFromSubtomo(subtomo)
            magn = subtomo.getAcquisition().getMagnification()
            rlnMicrographName = subtomo.getVolName()
            rlnCoordinateX = subtomo.getCoordinate3D().getX(BOTTOM_LEFT_CORNER)
            rlnCoordinateY = subtomo.getCoordinate3D().getY(BOTTOM_LEFT_CORNER)
            rlnCoordinateZ = subtomo.getCoordinate3D().getZ(BOTTOM_LEFT_CORNER)
            rlnImageName = subtomo.getFileName().replace(':' + MRC, '')
            rlnCtfImage = abspath(self._getCTFFileFromSubtomo(subtomo))
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

    @ staticmethod
    def _createStarTomoTable(isPyseg):
        pass
        # cols = RELION_TOMO_LABELS
        # # Pyseg post-rec only works if the magnification, pixel size and the prior angles aren't
        # # present in the star file
        # if isPyseg:
        #     cols = [TOMO_NAME,
        #             COORD_X,
        #             COORD_Y,
        #             COORD_Z,
        #             SUBTOMO_NAME,
        #             CTF_MISSING_WEDGE,
        #             ROT,
        #             TILT,
        #             PSI,
        #             SHIFTX,
        #             SHIFTY,
        #             SHIFTZ]
        # return Table(columns=cols)

    @ staticmethod
    def _getCTFFileFromSubtomo(subtomo):
        try:
            return subtomo.getCoordinate3D()._3dcftMrcFile.get()
        except:
            return 'Unavailable'

    @staticmethod
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
#############################################################################################################
    @staticmethod
    def _getTomogramStarFileLabels():
        return [
            TOMO_NAME,
            TILT_SERIES_NAME,
            CTFPLOTTER_FILE,
            IMOD_DIR,
            FRACTIONAL_DOSE,
            ACQ_ORDER_FILE,
            CULLED_FILE
        ]
    @staticmethod
    def _getSubtomogramStarFileLabels():
        return [
            TOMO_NAME,
            CTF_IMAGE,
            SUBTOMO_NAME,
            COORD_X,
            COORD_Y,
            COORD_Z,
            SHIFTX,
            SHIFTY,
            SHIFTZ,
            ROT,
            TILT,
            PSI,
            TILT_PRIOR,
            PSI_PRIOR,
            CLASS_NUMBER
        ]

    @staticmethod
    def _getCtfPlotterFile(tsId, ctfPlotterDir):
        return join(ctfPlotterDir, tsId, tsId + '.defocus')
