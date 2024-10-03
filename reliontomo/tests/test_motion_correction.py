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
from pyworkflow.tests import setupTestProject, DataSet
from pyworkflow.utils import magentaStr, yellowStr
from reliontomo import Plugin
from reliontomo.protocols import ProtRelionTomoMotionCorr
from tomo.objects import TomoAcquisition
from tomo.protocols import ProtImportTsMovies
from tomo.tests.test_base_centralized_layer import TestBaseCentralizedLayer

IS_RE_50 = Plugin.isRe50()


class TestRelion5MotionCorr(TestBaseCentralizedLayer):
    nTiltSeries = 2
    nTiltImgsPerTs = 3
    unbinnedTsDims = (7420, 7676)
    unbinnedSRate = 0.675
    # ACQUISITION PARAMS
    voltage = 300
    magnification = 105000
    sphericalAberration = 2.7
    amplitudeContrast = 0.1
    initialDose = 0
    dosePerFrame = 3.0
    tiltAxisAngle = 84.1
    # Only 3 images per TS, corresponding to the angles -60, 0, and 60 degrees with acq orders, respectively, 0, 39, 40
    accumDose = 120  # max acq order * dose per frame
    step = 60
    minAngle = -60
    maxAngle = 60

    @classmethod
    def setUpClass(cls):
        if IS_RE_50:
            setupTestProject(cls)
            cls.ds = DataSet.getDataSet('tomo-em')
            cls.testAcqObj = TomoAcquisition(
                voltage=cls.voltage,
                magnification=cls.magnification,
                sphericalAberration=cls.sphericalAberration,
                amplitudeContrast=cls.amplitudeContrast,
                samplingRate=cls.unbinnedSRate,
                doseInitial=cls.initialDose,
                dosePerFrame=cls.dosePerFrame,
                tiltAxisAngle=cls.tiltAxisAngle,
                angleMax=cls.maxAngle,
                angleMin=cls.minAngle,
                step=cls.step,
                accumDose=cls.accumDose
            )
            cls.importedTsM = cls._runImportTiltSeriesM()

    @classmethod
    def _runImportTiltSeriesM(cls, filesPattern='{TS}_{TO}_{TA}.mrc'):
        print(magentaStr("\n==> Importing data - tilt-series movies:"))
        protImport = cls.newProtocol(ProtImportTsMovies,
                                     filesPath=cls.ds.getFile('empiar'),
                                     filesPattern=filesPattern,
                                     voltage=cls.voltage,
                                     magnification=cls.magnification,
                                     sphericalAberration=cls.sphericalAberration,
                                     amplitudeContrast=cls.amplitudeContrast,
                                     samplingRate=cls.unbinnedSRate,
                                     doseInitial=cls.initialDose,
                                     dosePerFrame=cls.dosePerFrame,
                                     tiltAxisAngle=cls.tiltAxisAngle)
        cls.launchProtocol(protImport)
        return getattr(protImport, ProtImportTsMovies.OUTPUT_NAME, None)

    @classmethod
    def _runMotionCorr(cls, saveEvenOdd=False, binningFactor=1, patchX=1, patchY=1):
        print(magentaStr(f"\n==> Testing motioncorr:"
                         f"\n\t- Save odd/even = {saveEvenOdd}"
                         f"\n\t- Binning factor = {binningFactor}"
                         f"\n\t- PatchX, patchY = [{patchX}, {patchY}]"))
        protMotionCorr = cls.newProtocol(ProtRelionTomoMotionCorr,
                                         inputTiltSeriesM=cls.importedTsM,
                                         saveEvenOdd=saveEvenOdd,
                                         binningFactor=binningFactor,
                                         patchX=patchX,
                                         patchY=patchY)
        cls.launchProtocol(protMotionCorr)
        return protMotionCorr

    def _checkResults(self, protMotionCorr, saveEvenOdd=False, binningFactor=1):
        outTsSet = getattr(protMotionCorr, protMotionCorr._possibleOutputs.tiltSeries.name, None)
        outTsSetEven = getattr(protMotionCorr, protMotionCorr._possibleOutputs.tiltSeriesEven.name, None)
        outTsSetOdd = getattr(protMotionCorr, protMotionCorr._possibleOutputs.tiltSeriesOdd.name, None)
        expectedDims = [self.unbinnedTsDims[0] / binningFactor,
                        self.unbinnedTsDims[1] / binningFactor,
                        self.nTiltImgsPerTs]

        self._checkTiltSeries(outTsSet, binningFactor, expectedDims, saveEvenOdd)
        if saveEvenOdd:
            hasEvenOdd = False
            self._checkTiltSeries(outTsSetEven, binningFactor, expectedDims, hasEvenOdd)
            self._checkTiltSeries(outTsSetOdd, binningFactor, expectedDims, hasEvenOdd)
        else:
            self.assertIsNone(outTsSetEven)
            self.assertIsNone(outTsSetOdd)

    def _checkTiltSeries(self, tsSet, binningFactor, expectedDims, hasEvenOdd):
        self.checkTiltSeries(tsSet,
                             imported=True,  # There's no transformation matrix. Use this flag to skip that check
                             expectedSetSize=self.nTiltSeries,
                             expectedSRate=self.unbinnedSRate * binningFactor,
                             expectedDimensions=expectedDims,
                             testSetAcqObj=self.testAcqObj,
                             hasOddEven=hasEvenOdd,
                             anglesCountSet=self.nTiltImgsPerTs)

    def runTest(self, saveEvenOdd=False, binningFactor=1, patchX=1, patchY=1):
        if IS_RE_50:
            protMotionCorr = self._runMotionCorr(saveEvenOdd=saveEvenOdd,
                                                 binningFactor=binningFactor,
                                                 patchX=patchX,
                                                 patchY=patchY)
            self._checkResults(protMotionCorr,
                               saveEvenOdd=saveEvenOdd,
                               binningFactor=binningFactor)
        else:
            print(yellowStr('Relion 4 detected. Test for protocol "Motion correction" skipped.'))

    def testMotionCorr_01(self):
        self.runTest()

    def testMotionCorr_02(self):
        self.runTest(patchX=2, patchY=2)

    def testMotionCorr_03(self):
        self.runTest(saveEvenOdd=True, binningFactor=2)
