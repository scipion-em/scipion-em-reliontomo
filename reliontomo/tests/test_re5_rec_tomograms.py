# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors:     Scipion Team
# *
# * National Center of Biotechnology, CSIC, Spain
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
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************
import numpy as np
from imod.protocols import ProtImodImportTransformationMatrix
from imod.protocols.protocol_base import OUTPUT_TILTSERIES_NAME
from pyworkflow.tests import setupTestProject, DataSet
from pyworkflow.utils import magentaStr, yellowStr
from reliontomo import Plugin
from reliontomo.protocols import ProtRelion5TomoReconstruct
from reliontomo.protocols.protocol_re5_rec_tomogram import SINGLE_TOMO, ALL_TOMOS
from tomo.protocols import ProtImportTs, ProtImportTsCTF
from tomo.protocols.protocol_import_ctf import ImportChoice
from tomo.tests import RE4_STA_TUTO, DataSetRe4STATuto
from tomo.tests.test_base_centralized_layer import TestBaseCentralizedLayer

IS_RE_50 = Plugin.isRe50()


class TestRelionTomoRecTomograms(TestBaseCentralizedLayer):
    unbinnedSRate = DataSetRe4STATuto.unbinnedPixSize.value
    nTs = 2
    ts_03 = 'TS_03'
    ts_54 = 'TS_54'
    unbinnedWidth = 4000
    unbinnedHeight = 4000
    unbinnedThickness = 2000
    dims = [4000, 4000, 2000]

    @classmethod
    def setUpClass(cls):
        if IS_RE_50:
            setupTestProject(cls)
            cls.ds = DataSet.getDataSet(RE4_STA_TUTO)
            cls._runPreviousProtocols()

    @classmethod
    def _runPreviousProtocols(cls):
        cls.importedTs = cls._runImportTs()
        cls.importedCtfs = cls._runImportCtf()
        cls.tsWithAlignment = cls._runImportTrMatrix()

    @classmethod
    def _runImportTs(cls):
        print(magentaStr("\n==> Importing the tilt series:"))
        protImportTs = cls.newProtocol(ProtImportTs,
                                       filesPath=cls.ds.getFile(DataSetRe4STATuto.tsPath.value),
                                       filesPattern=DataSetRe4STATuto.tsPattern.value,
                                       exclusionWords=DataSetRe4STATuto.exclusionWordsTs03ts54.value,
                                       anglesFrom=2,  # From tlt file
                                       voltage=DataSetRe4STATuto.voltage.value,
                                       magnification=DataSetRe4STATuto.magnification.value,
                                       sphericalAberration=DataSetRe4STATuto.sphericalAb.value,
                                       amplitudeContrast=DataSetRe4STATuto.amplitudeContrast.value,
                                       samplingRate=cls.unbinnedSRate,
                                       doseInitial=DataSetRe4STATuto.initialDose.value,
                                       dosePerFrame=DataSetRe4STATuto.dosePerTiltImg.value,
                                       tiltAxisAngle=DataSetRe4STATuto.tiltAxisAngle.value)

        cls.launchProtocol(protImportTs)
        tsImported = getattr(protImportTs, 'outputTiltSeries', None)
        return tsImported

    @classmethod
    def _runImportCtf(cls):
        print(magentaStr("\n==> Importing the CTFs:"))
        protImportCtf = cls.newProtocol(ProtImportTsCTF,
                                        filesPath=cls.ds.getFile(DataSetRe4STATuto.tsPath.value),
                                        filesPattern=DataSetRe4STATuto.ctfPattern.value,
                                        importFrom=ImportChoice.CTFFIND.value,
                                        inputSetOfTiltSeries=cls.importedTs)
        cls.launchProtocol(protImportCtf)
        outputMask = getattr(protImportCtf, protImportCtf._possibleOutputs.CTFs.name, None)
        return outputMask

    @classmethod
    def _runImportTrMatrix(cls):
        print(magentaStr("\n==> Importing the TS' transformation matrices with IMOD:"))
        protImportTrMatrix = cls.newProtocol(ProtImodImportTransformationMatrix,
                                             filesPath=cls.ds.getFile(DataSetRe4STATuto.tsPath.value),
                                             filesPattern=DataSetRe4STATuto.transformPattern.value,
                                             inputSetOfTiltSeries=cls.importedTs)
        cls.launchProtocol(protImportTrMatrix)
        outTsSet = getattr(protImportTrMatrix, OUTPUT_TILTSERIES_NAME, None)
        return outTsSet

    def runRecTomograms(self, binnedPixSize=None, recTomoMode=None, tomoId=None, correctCtf=False):
        recTomoMsg = f'Reconstruct single tomo - tomoId = {tomoId}' if recTomoMode == SINGLE_TOMO \
            else "Reconstruct all the tomograms"
        print(magentaStr(f"\n==> {recTomoMsg}:"
                         f"\n\t- BinnedPixSize = {binnedPixSize:.3f} Å/px"
                         f"\n\t- CTF correction = {correctCtf}"))
        protRecTomos = self.newProtocol(ProtRelion5TomoReconstruct,
                                        inTsSet=self.tsWithAlignment,
                                        inCtfSet=self.importedCtfs,
                                        unbinnedWidth=self.unbinnedWidth,
                                        unbinnedHeight=self.unbinnedHeight,
                                        unbinnedThickness=self.unbinnedThickness,
                                        binnedPixSize=binnedPixSize,
                                        doCtfCorrection=correctCtf,
                                        recTomoMode=recTomoMode,
                                        tomoId=tomoId,
                                        numberOfThreads=1,
                                        binThreads=6)
        self.launchProtocol(protRecTomos)
        return getattr(protRecTomos, protRecTomos._possibleOutputs.tomograms.name, None)

    def _runTest(self, binningFactor=None, binnedPixSize=None, recTomoMode=None, tomoId=None, correctCtf=False):
        if IS_RE_50:
            # Run the protocol
            outTomos = self.runRecTomograms(binnedPixSize=binnedPixSize,
                                            recTomoMode=recTomoMode,
                                            tomoId=tomoId,
                                            correctCtf=correctCtf)
            # Check the results
            expectedDims = np.array(self.dims) / binningFactor
            self.checkTomograms(outTomos,
                                expectedSetSize=1 if recTomoMode == SINGLE_TOMO else self.nTs,
                                expectedSRate=binnedPixSize,
                                expectedDimensions=expectedDims.tolist(),
                                ctfCorrected=correctCtf)
        else:
            print(yellowStr('Relion 4 detected. Test for protocol "Reconstruct tomograms" skipped.'))

    def testRecTomos_01(self):
        binninFactor = 4
        binnedPixSize = self.unbinnedSRate * binninFactor
        recTomoMode = SINGLE_TOMO
        tomoId = self.ts_03
        correctCtf = True
        self._runTest(binningFactor=binninFactor,
                      binnedPixSize=binnedPixSize,
                      recTomoMode=recTomoMode,
                      tomoId=tomoId,
                      correctCtf=correctCtf)

    def testRecTomos_02(self):
        binninFactor = 10
        binnedPixSize = self.unbinnedSRate * binninFactor
        recTomoMode = SINGLE_TOMO
        tomoId = self.ts_54
        correctCtf = False
        self._runTest(binningFactor=binninFactor,
                      binnedPixSize=binnedPixSize,
                      recTomoMode=recTomoMode,
                      tomoId=tomoId,
                      correctCtf=correctCtf)

    def testRecTomos_03(self):
        binninFactor = 10
        binnedPixSize = self.unbinnedSRate * binninFactor
        recTomoMode = ALL_TOMOS
        correctCtf = True
        self._runTest(binningFactor=binninFactor,
                      binnedPixSize=binnedPixSize,
                      recTomoMode=recTomoMode,
                      correctCtf=correctCtf)
