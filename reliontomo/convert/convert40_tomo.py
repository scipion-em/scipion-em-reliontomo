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
from os import mkdir

from pwem import ALIGN_NONE, ALIGN_2D, ALIGN_PROJ
from pwem.emlib.image import ImageHandler
import pyworkflow.utils as pwutils
from pwem.objects import Transform
from pyworkflow.object import Float
from relion.convert.convert_base import WriterBase
from reliontomo.constants import TOMO_NAME, TILT_SERIES_NAME, CTFPLOTTER_FILE, IMOD_DIR, FRACTIONAL_DOSE, \
    ACQ_ORDER_FILE, CULLED_FILE, SUBTOMO_NAME, COORD_X, COORD_Y, COORD_Z, SHIFTX, SHIFTY, SHIFTZ, ROT, \
    TILT, PSI, TILT_PRIOR, PSI_PRIOR, CLASS_NUMBER, TOMO_PARTICLE_NAME, RANDOM_SUBSET, OPTICS_GROUP, SHIFTX_ANGST, \
    SHIFTY_ANGST, SHIFTZ_ANGST, CTF_IMAGE, FILE_NOT_FOUND
import pwem.convert.transformations as tfs
import numpy as np
from os.path import join
from pwem.convert.transformations import translation_from_matrix, euler_from_matrix
from relion.convert import Table, OpticsGroups
from reliontomo.objects import PseudoSubtomogram
from tomo.constants import BOTTOM_LEFT_CORNER
from tomo.objects import SubTomogram, Coordinate3D, TomoAcquisition


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
        precedentFileList = [tomo.getFileName() for tomo in precedentsSet]
        precedent = subtomo.getCoordinate3D().getVolume()
        # Get from the coordinate associated to the current subtomogram
        tsId = subtomo.getCoordinate3D().getTomoId()
        if not tsId:
            # Get from the tomogram in which that subtomo was picked
            if precedent:
                tsId = precedent.getTsId()
            # Match by filename
            if not tsId and subtomo.getVolName() in precedentFileList:
                tsId = precedentsSet[precedentFileList.index(subtomo.getVolName()) + 1].getTsId()

        if tsId:
            return tsId
        else:
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


class Reader:

    ALIGNMENT_LABELS = [
        "rlnOriginX",
        "rlnOriginY",
        "rlnOriginZ",
        "rlnAngleRot",
        "rlnAngleTilt",
        "rlnAnglePsi",
    ]

    def __init__(self, **kwargs):
        self.setParticleTransform = None
        self._shifts = None
        self._angles = None
        self._alignType = kwargs.get('alignType', ALIGN_NONE)
        self._pixelSize = kwargs.get('pixelSize', 1.0)

    def readPseudoSubtomgramsStarFile(self, starFile, precedents, outputSet):
        tomoTable = Table()
        tomoTable.read(starFile, tableName='particles')
        ih = ImageHandler()
        samplingRate = outputSet.getSamplingRate()
        precedentFileList, tsIdList = zip(*[[tomo.getFileName(), tomo.getTsId()] for tomo in precedents])
        precedentFileList = list(precedentFileList)
        tsIdList = list(tsIdList)
        opticsGroupStr = OpticsGroups.fromImages(outputSet).toString()

        for counter, row in enumerate(tomoTable):
            psubtomo = PseudoSubtomogram()
            coordinate3d = Coordinate3D()
            transform = Transform()
            origin = Transform()

            tomoName = row.get(TOMO_NAME, FILE_NOT_FOUND)
            tomoInd = tsIdList.index(tomoName)
            subtomoFilename = row.get(SUBTOMO_NAME, FILE_NOT_FOUND)
            coordinate3d.setVolume(precedents[tomoInd + 1])  # 3D coord must be referred to a volume to get its origin
            x = row.get(COORD_X, 0)
            y = row.get(COORD_Y, 0)
            z = row.get(COORD_Z, 0)
            tiltPrior = row.get(TILT_PRIOR, 0)
            psiPrior = row.get(PSI_PRIOR, 0)
            coordinate3d.setX(float(x), BOTTOM_LEFT_CORNER)
            coordinate3d.setY(float(y), BOTTOM_LEFT_CORNER)
            coordinate3d.setZ(float(z), BOTTOM_LEFT_CORNER)
            self.__setParticleTransformProj(psubtomo, row)

            psubtomo.setVolName(precedentFileList[tomoInd])
            psubtomo.setFileName(subtomoFilename)
            psubtomo.setCoordinate3D(coordinate3d)
            psubtomo.setTransform(transform)
            psubtomo.setAcquisition(TomoAcquisition())
            psubtomo.getAcquisition().opticsGroupInfo.set(opticsGroupStr)
            psubtomo.setClassId(row.get(CLASS_NUMBER, 0))
            psubtomo.setSamplingRate(samplingRate)
            psubtomo._tiltPriorAngle = Float(tiltPrior)
            psubtomo._psiPriorAngle = Float(psiPrior)

            # Set the origin and the dimensions of the current subtomogram
            x, y, z, n = ih.getDimensions(subtomoFilename)
            zDim, filename = self._manageIhDims(subtomoFilename, z, n)
            origin.setShifts(x / -2. * samplingRate, y / -2. * samplingRate, zDim / -2. * samplingRate)
            psubtomo.setOrigin(origin)

            # Set the pseudosubtomogram specific values
            psubtomo.setTomoParticleName(row.get(TOMO_PARTICLE_NAME, FILE_NOT_FOUND))
            psubtomo.setRandomSubset(row.get(RANDOM_SUBSET, 0))
            psubtomo.setOpticsGroupId(row.get(OPTICS_GROUP, 0))
            psubtomo.setShiftXAngst(row.get(SHIFTX_ANGST, 0))
            psubtomo.setShiftYAngst(row.get(SHIFTY_ANGST, 0))
            psubtomo.setShiftZAngst(row.get(SHIFTZ_ANGST, 0))
            psubtomo.setCtfImage(row.get(CTF_IMAGE, FILE_NOT_FOUND))

            # Add current subtomogram to the output set
            outputSet.append(psubtomo)

    @staticmethod
    def _manageIhDims(fileName, z, n):
        if fileName.endswith('.mrc') or fileName.endswith('.map'):
            fileName += ':mrc'
            if z == 1 and n != 1:
                zDim = n
            else:
                zDim = z
        else:
            zDim = z

        return zDim, fileName

    def setParticleTransform(self, particle, row):
        """ Set the transform values from the row. """

        if (self._alignType == ALIGN_NONE) or not row.hasAnyColumn(self.ALIGNMENT_LABELS):
            self.setParticleTransform = self.__setParticleTransformNone
        else:
            # Ensure the Transform object exists
            self._angles = np.zeros(3)
            self._shifts = np.zeros(3)

            particle.setTransform(Transform())

            if self._alignType == ALIGN_2D:
                self.setParticleTransform = self.__setParticleTransform2D
            elif self._alignType == ALIGN_PROJ:
                self.setParticleTransform = self.__setParticleTransformProj
            else:
                raise TypeError("Unexpected alignment type: %s"
                                % self._alignType)

        # Call again the modified function
        self.setParticleTransform(particle, row)

    @staticmethod
    def __setParticleTransformNone(particle, row):
        particle.setTransform(None)

    def __setParticleTransform2D(self, particle, row):
        angles = self._angles
        shifts = self._shifts

        shifts[0] = float(row.get(SHIFTX, 0))
        shifts[1] = float(row.get(SHIFTY, 0))
        angles[2] = float(row.get(PSI, 0))
        radAngles = -np.deg2rad(angles)
        M = tfs.euler_matrix(radAngles[0], radAngles[1], radAngles[2], 'szyz')
        M[:3, 3] = shifts[:3]
        particle.getTransform().setMatrix(M)

    def __setParticleTransformProj(self, particle, row):
        angles = self._angles
        shifts = self._shifts

        shifts[0] = float(row.get(SHIFTX, 0))
        shifts[1] = float(row.get(SHIFTY, 0))
        shifts[2] = float(row.get(SHIFTZ, 0))

        angles[0] = float(row.get(ROT, 0))
        angles[1] = float(row.get(TILT, 0))
        angles[2] = float(row.get(PSI, 0))

        radAngles = -np.deg2rad(angles)
        M = tfs.euler_matrix(radAngles[0], radAngles[1], radAngles[2], 'szyz')
        M[:3, 3] = -shifts[:3]
        M = np.linalg.inv(M)
        particle.getTransform().setMatrix(M)
