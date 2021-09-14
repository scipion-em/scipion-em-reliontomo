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
import csv
from pwem.emlib.image import ImageHandler
import pyworkflow.utils as pwutils
from relion.convert.convert_base import WriterBase
from reliontomo.constants import TOMO_NAME, TILT_SERIES_NAME, CTFPLOTTER_FILE, IMOD_DIR, FRACTIONAL_DOSE, \
    ACQ_ORDER_FILE, CULLED_FILE, SUBTOMO_NAME, COORD_X, COORD_Y, COORD_Z, SHIFTX, SHIFTY, SHIFTZ, ROT, \
    TILT, PSI, TILT_PRIOR, PSI_PRIOR, CLASS_NUMBER, TOMO_PARTICLE_NAME, RANDOM_SUBSET, OPTICS_GROUP, SHIFTX_ANGST, \
    SHIFTY_ANGST, SHIFTZ_ANGST, CTF_IMAGE
from scipion.install.funcs import mkdir
import numpy as np
from os.path import join
from pwem.convert.transformations import translation_from_matrix, euler_from_matrix
from relion.convert import Table
from tomo.constants import BOTTOM_LEFT_CORNER


class Writer(WriterBase):
    """ Helper class to convert from Scipion SetOfTomograms and SetOfSubTomograms to star files ."""

    def writeSetOfTomograms(self, tomoSet, outStarFileName, prot=None, tsSet=None,
                            ctfPlotterParentDir=None, eTomoParentDir=None):
        tomoTable = Table(columns=self._getTomogramStarFileLabels())
        tsIdFromTsSet = [ts.getTsId() for ts in tsSet]
        for tomo in tomoSet:
            tsId = tomo.getTsId()
            ts = tsSet[tsIdFromTsSet.index(tsId) + 1]
            tomoTable.addRow(tsId,                                                # _rlnTomoName #1
                             ts.getFirstItem().getFileName() + ':mrc',            # _rlnTomoTiltSeriesName #2
                             self._getCtfPlotterFile(tsId, ctfPlotterParentDir),  # _rlnTomoImportCtfPlotterFile #3
                             join(eTomoParentDir, tsId),                          # _rlnTomoImportImodDir #4
                             ts.getAcquisition().getDosePerFrame(),               # _rlnTomoImportFractionalDose #5
                             self._genOrderListFile(prot, ts),                    # _rlnTomoImportOrderList #6
                             self._genCulledFileName(prot, tsId)                  # _rlnTomoImportCulledFile #7
                             )

        # Write the STAR file
        tomoTable.write(outStarFileName)

    def writeSetOfSubtomograms(self, subtomoSet, subtomosStar):
        currentTomo = ''
        MRC = 'mrc'
        ih = ImageHandler()
        tomoTable = Table(columns=self._getSubtomogramStarFileLabels())
        tmpDir = pwutils.getParentFolder(subtomosStar)
        precedentsSet = subtomoSet.getCoordinates3D().get().getPrecedents()
        for subtomo in subtomoSet:
            if pwutils.getExt(subtomo.getFileName().replace(':' + MRC, '')) != '.' + MRC:
                mrcDir = join(tmpDir, pwutils.removeBaseExt(subtomo.getVolName()))
                if currentTomo != subtomo.getVolName():
                    mkdir(mrcDir)
                mrcFile = join(mrcDir, pwutils.replaceBaseExt(subtomo.getFileName(), MRC))
                ih.convert(subtomo.getFileName(), mrcFile)
            angles, shifts = self._getTransformInfoFromSubtomo(subtomo)
            rlnAngleTilt = angles[1]
            rlnAnglePsi = angles[2]
            # Add row to the table which will be used to generate the STAR file
            tomoTable.addRow(
                self._getPrecedentTsIdFromSubtomo(subtomo, precedentsSet),   # _rlnTomoName #1
                subtomo.getFileName().replace(':' + MRC, ''),                # _rlnImageName #2
                subtomo.getCoordinate3D().getX(BOTTOM_LEFT_CORNER),          # _rlnCoordinateX #3
                subtomo.getCoordinate3D().getY(BOTTOM_LEFT_CORNER),          # _rlnCoordinateY #4
                subtomo.getCoordinate3D().getZ(BOTTOM_LEFT_CORNER),          # _rlnCoordinateZ #5
                shifts[0],                                                   # _rlnOriginX #6
                shifts[1],                                                   # _rlnOriginY #7
                shifts[2],                                                   # _rlnOriginZ #8
                angles[0],                                                   # _rlnAngleRot #9
                rlnAngleTilt,                                                # _rlnAngleTilt #10
                rlnAnglePsi,                                                 # _rlnAnglePsi #11
                getattr(subtomo, '_tiltPriorAngle', rlnAngleTilt),           # _rlnAngleTiltPrior #12
                getattr(subtomo, '_psiPriorAngle', rlnAnglePsi),             # _rlnAnglePsiPrior #13
                getattr(subtomo, '_getClassId', 0)                           # _rlnClassNumber #14
            )

        # Write the STAR file
        tomoTable.write(subtomosStar)

    @staticmethod
    def _getPrecedentTsIdFromSubtomo(subtomo, precedentsSet):
        for tomo in precedentsSet:
            if subtomo.getVolName() in tomo.getFileName():
                return tomo.getTsId()

        raise Exception('Not able to read TsId for subtomogram %s' % subtomo.getVolName())

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
    def _getPseudoSubtomogramStarFileLabels():
        return [
            TOMO_NAME,
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
            CLASS_NUMBER,
            TOMO_PARTICLE_NAME,
            RANDOM_SUBSET,
            OPTICS_GROUP,
            SHIFTX_ANGST,
            SHIFTY_ANGST,
            SHIFTZ_ANGST,
            CTF_IMAGE
        ]

    @staticmethod
    def _getCtfPlotterFile(tsId, ctfPlotterDir):
        return join(ctfPlotterDir, tsId, tsId + '.defocus')

    @staticmethod
    def _genOrderListFile(prot, ts):
        """The order file expected by Relion is A 2-column, comma-separated file with the frame-order list
        of the tilt series, where the first column is the frame (image) number (starting at 1) and the second
        column is the tilt angle (in degrees).
        :param prot: current protocol object
        :param ts: TiltSeries object"""
        outputFilename = prot._getExtraPath(ts.getTsId() + '_order_list.csv')
        tiList = [ti.clone() for ti in ts]
        ind = np.argsort([ti.getTiltAngle() for ti in tiList])  # Indices to get the data sorted by acqOrder
        with open(outputFilename, mode='w') as acqOrderFile:
            acqOrderFileWriter = csv.writer(acqOrderFile, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
            acqOrderList = [ti.getAcquisitionOrder() for ti in tiList]
            if min(acqOrderList) == 0:
                [acqOrderFileWriter.writerow([tiList[i].getAcquisitionOrder() + 1, tiList[i].getTiltAngle()]) for i in ind]
            else:
                [acqOrderFileWriter.writerow([tiList[i].getAcquisitionOrder(), tiList[i].getTiltAngle()]) for i in ind]

        return outputFilename

    @staticmethod
    def _genCulledFileName(prot, tsId):
        return prot._getExtraPath(tsId + '_culled.mrc:mrc')
