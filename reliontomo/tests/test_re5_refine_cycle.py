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
from os.path import exists
from typing import Type
from imod.protocols import ProtImodImportTransformationMatrix
from imod.protocols.protocol_base import OUTPUT_TILTSERIES_NAME
from pwem.protocols import ProtImportMask, ProtImportVolumes
from pwem.protocols.protocol_import.masks import ImportMaskOutput
from pyworkflow.tests import setupTestProject, DataSet
from pyworkflow.utils import magentaStr, yellowStr
from reliontomo import Plugin
from reliontomo.constants import OUT_PARTICLES_STAR, IN_TOMOS_STAR, FILE_NOT_FOUND, POSTPROCESS_STAR_FIELD
from reliontomo.convert.convertBase import getTransformInfoFromCoordOrSubtomo
from reliontomo.protocols import ProtImportCoordinates3DFromStar, ProtRelion5ExtractSubtomos, \
    ProtRelion5ReconstructParticle, ProtRelionDeNovoInitialModel, ProtRelionRefineSubtomograms, \
    ProtRelion3DClassifySubtomograms, ProtRelionPostProcess
from reliontomo.protocols.protocol_edit_particles_star import OPERATION_LABELS, LABELS_TO_OPERATE_WITH, \
    ProtRelionEditParticlesStar
from reliontomo.tests import DataSetRe4Tomo
from reliontomo.utils import genEnumParamDict
from tomo.constants import TR_SCIPION
from tomo.objects import SetOfCoordinates3D
from tomo.protocols import ProtImportTs, ProtImportTsCTF, ProtImportTomograms
from tomo.protocols.protocol_import_ctf import ImportChoice
from tomo.protocols.protocol_import_tomograms import OUTPUT_NAME
from tomo.tests import RE4_STA_TUTO, DataSetRe4STATuto
from tomo.tests.test_base_centralized_layer import TestBaseCentralizedLayer


IS_RE_50 = Plugin.isRe50()


class TestRelion5RefineCycleBase(TestBaseCentralizedLayer):
    ds = None
    tsWithAlignment = None
    importedTs = None
    importedCtfs = None
    importedTomos = None
    importedCoords = None
    extractedParticles = None
    nParticles = DataSetRe4STATuto.nCoordsTotal.value
    nImgs = 41
    unbinnedSRate = DataSetRe4STATuto.unbinnedPixSize.value
    croppedBoxSizeBin6 = DataSetRe4STATuto.croppedBoxSizeBin4.value
    boxSizeBin6 = DataSetRe4STATuto.boxSizeBin4.value
    croppedBoxSizeBin2 = DataSetRe4STATuto.croppedBoxSizeBin2.value
    boxSizeBin2 = DataSetRe4STATuto.boxSizeBin2.value
    tsIds = ['TS_03', 'TS_54']
    binFactor4 = 4
    binFactor6 = 6
    binFactor2 = 2
    protExtract = None
    importedRef = None
    inReParticles = None
    symmetry = DataSetRe4STATuto.symmetry.value

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
        cls.importedTomos = cls._runImportTomograms()
        cls.importedCoords = cls._runImportCoordinatesFromStar()

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
    def _runImportFscMask(cls):
        print(magentaStr("\n==> Importing the FSC mask:"))
        protImportMask = cls.newProtocol(ProtImportMask,
                                         maskPath=cls.ds.getFile(DataSetRe4Tomo.maskFscBin2.name),
                                         samplingRate=cls.unbinnedSRate)
        protImportMask.setObjLabel('Import FSC mask, bin 2')
        cls.launchProtocol(protImportMask)
        outputMask = getattr(protImportMask, ImportMaskOutput.outputMask.name, None)
        return outputMask

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

    @classmethod
    def _runImportTomograms(cls):
        print(magentaStr("\n==> Importing the tomograms:"))
        protImportTomos = cls.newProtocol(ProtImportTomograms,
                                          filesPath=cls.ds.getFile(DataSetRe4STATuto.tsPath.value),
                                          filesPattern=DataSetRe4STATuto.tomosPattern.value,
                                          samplingRate=DataSetRe4STATuto.sRateBin4.value)  # Bin 4
        cls.launchProtocol(protImportTomos)
        outTomos = getattr(protImportTomos, OUTPUT_NAME, None)
        return outTomos

    @classmethod
    def _runImportCoordinatesFromStar(cls):
        print(magentaStr("\n==> Importing the coordinates with Relion:"))
        protImportCoords = cls.newProtocol(ProtImportCoordinates3DFromStar,
                                           starFile=cls.ds.getFile(DataSetRe4STATuto.coordsStarSubset.value),
                                           inTomos=cls.importedTomos,
                                           samplingRate=cls.unbinnedSRate)
        cls.launchProtocol(protImportCoords)
        outCoords = getattr(protImportCoords, ProtImportCoordinates3DFromStar._possibleOutputs.coordinates.name, None)
        return outCoords

    @classmethod
    def _runExtractSubtomos(cls, inParticles=None, inputCtfTs=None, inputTS=None, binningFactor=None,
                            boxSize=None, croppedBoxSize=None, gen2dParticles=False):
        extractMsg = 'Extract' if type(inParticles) is SetOfCoordinates3D else 'Re-extract'
        partTypeMsg = '2D' if gen2dParticles else '3D'
        print(magentaStr(f"\n==> {extractMsg}ing the particles from the tilt-series:\n"
                         f"\t- Particles type selected: {partTypeMsg}"))
        protExtract = cls.newProtocol(ProtRelion5ExtractSubtomos,
                                      inputCtfTs=inputCtfTs,
                                      inReParticles=inParticles,
                                      inputTS=inputTS,
                                      binningFactor=binningFactor,
                                      boxSize=boxSize,
                                      croppedBoxSize=croppedBoxSize,
                                      write2dStacks=gen2dParticles,
                                      numberOfMpi=2,  # There are 2 tomograms in the current dataset
                                      binThreads=4)
        protExtract.setObjLabel(f"{extractMsg} subtomos {partTypeMsg}")
        return cls.launchProtocol(protExtract)

    @classmethod
    def _runImportReference(cls, filePath=None, samplingRate=None):
        print(magentaStr("\n==> Importing the reference volume:"))
        protImportRef = cls.newProtocol(ProtImportVolumes,
                                        filesPath=filePath,
                                        samplingRate=samplingRate)
        cls.launchProtocol(protImportRef)
        return getattr(protImportRef, ProtImportVolumes._possibleOutputs.outputVolume.name, None)

    def _checkRe4Metadata(self, pSubtomoSet, noParticles=-1, tomogramsFile=None, particlesFile=None,
                          trajectoriesFile=None,
                          referenceFscFile=None, relionBinning=None, are2dStacks=-1, tsSamplingRate=None):
        self.assertEqual(noParticles, pSubtomoSet.getNReParticles())
        self.assertEqual(self.unbinnedSRate, pSubtomoSet.getTsSamplingRate())
        self.assertEqual(tomogramsFile, pSubtomoSet.getTomogramsStar())
        self.assertTrue(particlesFile, pSubtomoSet.getParticlesStar())
        if trajectoriesFile:
            self.assertEqual(trajectoriesFile, pSubtomoSet.getTrajectoriesStar())
        else:
            # It may exist in the optimisation set, but as an empty string
            self.assertTrue(pSubtomoSet.getTrajectoriesStar() in ["", None])
        if referenceFscFile:
            self.assertEqual(referenceFscFile, pSubtomoSet.getReferenceFsc())
        else:
            # It may exist in the optimisation set, but as an empty string
            self.assertTrue(pSubtomoSet.getReferenceFsc() in ["", None])
        self.assertEqual(pSubtomoSet.getTsSamplingRate(), tsSamplingRate)
        self.assertEqual(pSubtomoSet.getCurrentSamplingRate(), tsSamplingRate * relionBinning)
        self.assertEqual(pSubtomoSet.getRelionBinning(), relionBinning)
        self.assertEqual(pSubtomoSet.are2dStacks(), are2dStacks)
        self.assertTrue(pSubtomoSet.areRe5Particles())

    def _checkPseudosubtomograms(self, inCoords, outSubtomos, expectedSetSize=-1, expectedSRate=-1, expectedBoxSize=-1,
                                 convention=TR_SCIPION, orientedParticles=False, are2dStacks=None):
        # This way the box size check will be skipped in checkExtractedSubtomos and the corresponding box size will
        # be checked below here
        expectedExtension = '.mrcs' if are2dStacks else '.mrc'
        # Check the SubTomogram part, as Relionparticles inherit from them
        self.checkExtractedSubtomos(inCoords, outSubtomos,
                                    expectedSetSize=expectedSetSize,
                                    expectedSRate=expectedSRate,
                                    expectedBoxSize=expectedBoxSize,
                                    convention=convention,
                                    orientedParticles=orientedParticles,
                                    expectedExtension=expectedExtension,
                                    isStack2d=are2dStacks,
                                    noImgs=self.nImgs)

        # Check some remaining specific relionParticle attributes
        if are2dStacks:
            self.assertTrue(outSubtomos.are2dStacks())
            for pSubtomo in outSubtomos:
                self.assertGreater(len(pSubtomo.getVisibleFrames()), 1)  # List of 0s and 1os if filled, 0 if not
                self.assertEqual(pSubtomo.getCtfFile(), FILE_NOT_FOUND)  # Only generated for 3d particles
                self.assertTrue(pSubtomo.getTsId() in self.tsIds)
        else:
            self.assertFalse(outSubtomos.are2dStacks())
            for pSubtomo in outSubtomos:
                self.assertEqual(len(pSubtomo.getVisibleFrames()), 1)
                self.assertTrue(exists(pSubtomo.getCtfFile()))
                self.assertTrue(pSubtomo.getTsId() in self.tsIds)


class TestRelion5TomoExtractSubtomos(TestRelion5RefineCycleBase):

    def _runTestExtractSubtomos(self, inParticles=None, inputCtfTs=None, inputTS=None, binningFactor=None,
                                boxSize=None, croppedBoxSize=None, are2dParticles=False, isReExtraction=False):
        protExtract = self._runExtractSubtomos(inParticles=inParticles,
                                               inputCtfTs=inputCtfTs,
                                               inputTS=inputTS,
                                               binningFactor=binningFactor,
                                               boxSize=boxSize,
                                               croppedBoxSize=croppedBoxSize,
                                               gen2dParticles=are2dParticles)
        outParticles = getattr(protExtract, ProtRelion5ExtractSubtomos._possibleOutputs.relionParticles.name, None)
        # Check RelionTomoMetadata: both particles and tomograms files are generated
        self._checkRe4Metadata(outParticles,
                               noParticles=self.nParticles,
                               tomogramsFile=protExtract._getExtraPath(IN_TOMOS_STAR),
                               particlesFile=protExtract._getExtraPath(OUT_PARTICLES_STAR),
                               trajectoriesFile=None,
                               referenceFscFile=None,
                               relionBinning=binningFactor,
                               are2dStacks=are2dParticles,
                               tsSamplingRate=self.unbinnedSRate)
        # Check that the projected coordinates have been generated
        fiducials = getattr(protExtract, ProtRelion5ExtractSubtomos._possibleOutputs.projected2DCoordinates.name, None)
        if isReExtraction:
            self.assertIsNone(fiducials, msg='The projected coordinates wer not expected to be genrated.')
        else:
            self.assertIsNotNone(fiducials, msg='The projected coordinates were not generated.')
        # Check the pseudo-subtomograms
        unbinnedPixSize = self.unbinnedSRate
        nParticles = self.nParticles
        self._checkPseudosubtomograms(self.importedCoords, outParticles,
                                      expectedSetSize=nParticles,
                                      expectedSRate=unbinnedPixSize * binningFactor,
                                      expectedBoxSize=croppedBoxSize,
                                      orientedParticles=True,  # The imported coordinates are oriented
                                      are2dStacks=are2dParticles)
        return outParticles

    def testExtractSubtomos3d(self):
        if IS_RE_50:
            self._runTestExtractSubtomos(inParticles=self.importedCoords,
                                         inputCtfTs=self.importedCtfs,
                                         inputTS=self.tsWithAlignment,
                                         binningFactor=self.binFactor6,
                                         boxSize=self.boxSizeBin6,
                                         croppedBoxSize=self.croppedBoxSizeBin6,
                                         are2dParticles=False)
        else:
            print(yellowStr('Relion 4 detected. Test for protocol "Extract subtomos" skipped.'))

    def testExtractSubtomos2d(self):
        if IS_RE_50:
            inputCtfTs = self.importedCtfs
            inputTS = self.tsWithAlignment
            extractedParticles = self._runTestExtractSubtomos(inParticles=self.importedCoords,
                                                              inputCtfTs=inputCtfTs,
                                                              inputTS=inputTS,
                                                              binningFactor=self.binFactor6,
                                                              boxSize=self.boxSizeBin6,
                                                              croppedBoxSize=self.croppedBoxSizeBin6,
                                                              are2dParticles=True)
            # Run a re-extraction
            self._runTestExtractSubtomos(inParticles=extractedParticles,
                                         inputCtfTs=inputCtfTs,
                                         inputTS=inputTS,
                                         binningFactor=self.binFactor2,
                                         boxSize=self.boxSizeBin2,
                                         croppedBoxSize=self.croppedBoxSizeBin2,
                                         are2dParticles=True,
                                         isReExtraction=True)

            # Run a re-extraction without input TS nor CTFs
            self._runTestExtractSubtomos(inParticles=extractedParticles,
                                         binningFactor=self.binFactor2,
                                         boxSize=self.boxSizeBin2,
                                         croppedBoxSize=self.croppedBoxSizeBin2,
                                         are2dParticles=True,
                                         isReExtraction=True)
        else:
            print(yellowStr('Relion 4 detected. Test for protocol "Extract subtomos" skipped.'))


class TestRelion5TomoRecParticleFromTs(TestRelion5RefineCycleBase):

    def testRecParticleFromTS_01(self):
        self._runTest(are2dStacks=True, hasHalves=False)  # Halves are generated after having run a refinement

    def testRecParticleFromTS_02(self):
        self._runTest(are2dStacks=False, hasHalves=False)

    def _runTest(self, are2dStacks=True, hasHalves=None):
        if IS_RE_50:
            binningFactor = self.binFactor6
            bozSize = self.boxSizeBin6
            croppedBoxSize = self.croppedBoxSizeBin6
            protExtract = self._runExtractSubtomos(inParticles=self.importedCoords,
                                                   inputCtfTs=self.importedCtfs,
                                                   inputTS=self.tsWithAlignment,
                                                   binningFactor=binningFactor,
                                                   boxSize=bozSize,
                                                   croppedBoxSize=croppedBoxSize,
                                                   gen2dParticles=are2dStacks)
            extractedPars = getattr(protExtract, ProtRelion5ExtractSubtomos._possibleOutputs.relionParticles.name, None)
            protRecPartFromTS = self._recParticleFromTS(extractedPars,
                                                        binningFactor=binningFactor,
                                                        boxSize=bozSize,
                                                        croppedBoxSize=croppedBoxSize,
                                                        are2dStacks=are2dStacks)
            # Check the results
            self._checkResults(protRecPartFromTS,
                               binningFactor=binningFactor,
                               expectedBoxSize=croppedBoxSize,
                               hasHalves=hasHalves)
        else:
            print(yellowStr('Relion 4 detected. Test for protocol "Reconstruct particle from TS" skipped.'))

    @classmethod
    def _recParticleFromTS(cls, extractedParticles, binningFactor=None, boxSize=None, croppedBoxSize=None,
                           are2dStacks=True):
        partTypeStr = '2D' if are2dStacks else '3D'
        print(magentaStr(f"\n==> Reconstructing the particle from the TS:"
                         f"\n\t- Using {partTypeStr} stacks."
                         f"\n\t- Binning factor = {binningFactor}."))
        protRecPartFromTS = cls.newProtocol(ProtRelion5ReconstructParticle,
                                            inReParticles=extractedParticles,
                                            binningFactor=binningFactor,
                                            boxSize=boxSize,
                                            croppedBoxSize=croppedBoxSize,
                                            symmetry='C6',  # HIV symmetry
                                            numberOfMpi=2,  # 2 TS in the dataset used
                                            binThreads=4)
        protRecPartFromTS.setObjLabel(f'Rec particle from TS, {partTypeStr} stacks')
        cls.launchProtocol(protRecPartFromTS)
        return protRecPartFromTS

    def _checkResults(self, protRecPartFromTS, binningFactor=None, expectedBoxSize=None, hasHalves=None):
        # Check the reconstructed volume
        recVol = getattr(protRecPartFromTS, ProtRelion5ReconstructParticle._possibleOutputs.average.name, None)
        self.checkAverage(recVol,
                          expectedSRate=self.unbinnedSRate * binningFactor,
                          expectedBoxSize=expectedBoxSize,
                          hasHalves=hasHalves)


class TestRelion5TomoEditStar(TestRelion5RefineCycleBase):
    editStarOperationDict = genEnumParamDict(OPERATION_LABELS)
    editStarLabelsDict = genEnumParamDict(LABELS_TO_OPERATE_WITH)
    editTestsTol = 0.01

    @classmethod
    def _runPreviousProtocols(cls):
        super()._runPreviousProtocols()
        cls.protExtract = cls._runExtractSubtomos(inParticles=cls.importedCoords,
                                                  inputCtfTs=cls.importedCtfs,
                                                  inputTS=cls.tsWithAlignment,
                                                  binningFactor=cls.binFactor6,
                                                  boxSize=cls.boxSizeBin6,
                                                  croppedBoxSize=cls.croppedBoxSizeBin6,
                                                  gen2dParticles=True)
        cls.inReParticles = getattr(cls.protExtract, cls.protExtract._possibleOutputs.relionParticles.name, None)

    def testEditStar_shiftCenter(self):
        if IS_RE_50:
            print(magentaStr("\n==> Performing re-centering of particles:"))
            shiftX_pix = 4
            shiftY_pix = 2
            shiftZ_pix = 3
            protEditStar = self.newProtocol(ProtRelionEditParticlesStar,
                                            inReParticles=self.inReParticles,
                                            doRecenter=True,
                                            shiftX=shiftX_pix,
                                            shiftY=shiftY_pix,
                                            shiftZ=shiftZ_pix)
            protEditStar.setObjLabel('Particles re-center')
            self.launchProtocol(protEditStar)
            # Check the results
            inPSubtomos = protEditStar.getInputParticles()
            outPSubtomos = getattr(protEditStar, protEditStar._possibleOutputs.relionParticles.name, None)
            sRate = self.unbinnedSRate * self.binFactor6
            for inPSubtomo, outPSubtomo in zip(inPSubtomos, outPSubtomos):
                isx, isy, isz = self._getShiftsFromPSubtomogram(inPSubtomo)
                osx, osy, osz = self._getShiftsFromPSubtomogram(outPSubtomo)
                self.assertTrue(abs((isx - shiftX_pix * sRate) - osx) < self.editTestsTol)
                self.assertTrue(abs((isy - shiftY_pix * sRate) - osy) < self.editTestsTol)
                self.assertTrue(abs((isz - shiftZ_pix * sRate) - osz) < self.editTestsTol)
        else:
            print(yellowStr('Relion 4 detected. Test for protocol "Edit particles star" skipped.'))

    def testEditStar_removeDuplicates(self):
        if IS_RE_50:
            print(magentaStr("\n==> Performing particle removal:"))
            protEditStar = self.newProtocol(ProtRelionEditParticlesStar,
                                            inReParticles=self.inReParticles,
                                            doRemoveDuplicates=True,
                                            minDistPartRemoval=70  # In angstroms
                                            )
            protEditStar.setObjLabel('Particle removal')
            self.launchProtocol(protEditStar)
            # Check the results
            expectedSetSize = 176
            outParticles = getattr(protEditStar, protEditStar._possibleOutputs.relionParticles.name, None)
            binningFactor = self.binFactor6
            unbinnedPixSize = self.unbinnedSRate
            self._checkRe4Metadata(outParticles,
                                   noParticles=expectedSetSize,
                                   particlesFile=protEditStar._getExtraPath(OUT_PARTICLES_STAR),
                                   tomogramsFile=protEditStar.getInputParticles().getTomogramsStar(),
                                   relionBinning=binningFactor,
                                   are2dStacks=True,
                                   tsSamplingRate=unbinnedPixSize
                                   )
            self._checkPseudosubtomograms(self.importedCoords, outParticles,
                                          expectedSetSize=expectedSetSize,
                                          expectedSRate=unbinnedPixSize * binningFactor,
                                          expectedBoxSize=self.croppedBoxSizeBin6,
                                          orientedParticles=True,  # The imported coordinates are oriented
                                          are2dStacks=True)
        else:
            print(yellowStr('Relion 4 detected. Test for protocol "Edit particles star" skipped.'))

    @classmethod
    def _editStar_shiftCenter(cls) -> Type[ProtRelionEditParticlesStar]:
        print(magentaStr("\n==> Perform centering of particles:"))
        protEditStar = cls.newProtocol(ProtRelionEditParticlesStar,
                                       inReParticles=cls.inReParticles,
                                       doRecenter=True,
                                       shiftX=4,
                                       shiftY=2,
                                       shiftZ=3)
        protEditStar.setObjLabel('Particles re-center')
        cls.launchProtocol(protEditStar)
        return protEditStar

    @classmethod
    def _getShiftsFromPSubtomogram(cls, pSubtomo):
        _, shifts = getTransformInfoFromCoordOrSubtomo(pSubtomo, pSubtomo.getSamplingRate())
        return shifts[0], shifts[1], shifts[2]


class TestRelion5TomoGenInitialModel(TestRelion5RefineCycleBase):

    def testInitialModel(self):
        if IS_RE_50:
            recVol = self._genInitialModel()
            self.checkAverage(recVol,
                              expectedSRate=self.unbinnedSRate * self.binFactor6,
                              expectedBoxSize=self.croppedBoxSizeBin6,
                              hasHalves=False)
        else:
            print(yellowStr('Relion 4 detected. Test for protocol "Generate initial volume" skipped.'))

    @classmethod
    def _genInitialModel(cls, are2dStacks=True):
        protExtract = cls._runExtractSubtomos(
            inParticles=cls.importedCoords,
            inputCtfTs=cls.importedCtfs,
            inputTS=cls.tsWithAlignment,
            binningFactor=cls.binFactor6,
            boxSize=cls.boxSizeBin6,
            croppedBoxSize=cls.croppedBoxSizeBin6,
            gen2dParticles=are2dStacks)
        extractedPars = getattr(protExtract, ProtRelion5ExtractSubtomos._possibleOutputs.relionParticles.name, None)

        print(magentaStr("\n==> Generating the a de novo 3D initial model:"))
        protInitialModel = cls.newProtocol(ProtRelionDeNovoInitialModel,
                                           inReParticles=extractedPars,
                                           nVdamMiniBatches=20,
                                           maskDiameter=500,
                                           symmetry=cls.symmetry,
                                           priorWidthTiltAngle=10,
                                           doInC1AndApplySymLater=False,
                                           pooledSubtomos=3,
                                           doGpu=True,
                                           gpusToUse='0',
                                           numberOfMpi=1,
                                           binThreads=6)
        cls.launchProtocol(protInitialModel)
        return getattr(protInitialModel, ProtRelionDeNovoInitialModel._possibleOutputs.average.name, None)


class TestRelionTomo5Classify3d(TestRelion5RefineCycleBase):
    nClasses = 3

    @classmethod
    def _runPreviousProtocols(cls):
        super()._runPreviousProtocols()
        cls.importedRef = cls._runImportReference(filePath=cls.ds.getFile(DataSetRe4STATuto.recParticleBin6.name),
                                                  samplingRate=cls.unbinnedSRate * cls.binFactor6)
        cls.protExtract = cls._runExtractSubtomos(inParticles=cls.importedCoords,
                                                  inputCtfTs=cls.importedCtfs,
                                                  inputTS=cls.tsWithAlignment,
                                                  binningFactor=cls.binFactor6,
                                                  boxSize=cls.boxSizeBin6,
                                                  croppedBoxSize=cls.croppedBoxSizeBin6,
                                                  gen2dParticles=True)
        cls.extractedParts = getattr(cls.protExtract,
                                     ProtRelion5ExtractSubtomos._possibleOutputs.relionParticles.name,
                                     None)

    def testCl3d(self):
        if IS_RE_50:
            protCl3d = self._run3dClassify()
            self._checkCl3dResults(protCl3d)
        else:
            print(yellowStr('Relion 4 detected. Test for protocol "3d classification" skipped.'))

    def testCl3dWithAlignment(self):
        if IS_RE_50:
            protCl3d = self._run3dClassify(doAlingment=True)
            self._checkCl3dResults(protCl3d, onlyClassify=False)
        else:
            print(yellowStr('Relion 4 detected. Test for protocol "3d classification" skipped.'))

    @classmethod
    def _run3dClassify(cls, doAlingment=False):
        paramsDict = {'inReParticles': cls.extractedParts,
                      'referenceVolume': cls.importedRef,
                      'numberOfClasses': cls.nClasses,
                      'initialLowPassFilterA': 60,
                      'symmetry': cls.symmetry,
                      'maskDiameter': 500,
                      'pooledSubtomos': 6,
                      'numberOfMpi': 3,
                      'binThreads': 3}

        if doAlingment:
            paramsDict['doImageAlignment'] = True
            paramsDict['nIterations'] = 1  # Prevent error "No orientation was found as better than any other"
            # because of the test dataset
            paramsDict['doGpu'] = True
            paramsDict['gpusToUse'] = '0'
            label = 'cl3d with alignment'
            print(magentaStr("\n==> Classifying and aligning the pseudo-subtomograms:"))
        else:
            paramsDict['nIterations'] = 3
            label = 'cl3d'
            print(magentaStr("\n==> Classifying the pseudo-subtomograms:"))

        protCl3d = cls.newProtocol(ProtRelion3DClassifySubtomograms, **paramsDict)
        protCl3d.setObjLabel(label)
        cls.launchProtocol(protCl3d)
        return protCl3d

    def _checkCl3dResults(self, protCl3d, are2dStacks=True, onlyClassify=True):
        relionParticles = getattr(protCl3d, ProtRelion3DClassifySubtomograms._possibleOutputs.relionParticles.name,
                                  None)
        outClasses = getattr(protCl3d, ProtRelion3DClassifySubtomograms._possibleOutputs.classes.name, None)
        # Check RelionTomoMetadata: only the particles file is generated
        self._checkRe4Metadata(relionParticles,
                               noParticles=self.nParticles,
                               tomogramsFile=self.protExtract._getExtraPath(IN_TOMOS_STAR),
                               particlesFile=protCl3d._getExtraPath(OUT_PARTICLES_STAR),
                               trajectoriesFile=None,
                               referenceFscFile=None,
                               relionBinning=self.binFactor6,
                               are2dStacks=are2dStacks,
                               tsSamplingRate=self.unbinnedSRate)

        # Check the set of pseudosubtomograms
        expectedSetSize = self.nParticles
        expectedSRate = self.unbinnedSRate * self.binFactor6
        expectedBoxSize = self.croppedBoxSizeBin6
        if onlyClassify:
            self._checkPseudosubtomograms(self.importedCoords, relionParticles,
                                          expectedSetSize=expectedSetSize,
                                          expectedSRate=expectedSRate,
                                          expectedBoxSize=expectedBoxSize,
                                          orientedParticles=True,
                                          are2dStacks=are2dStacks)
        else:
            self.checkRefinedSubtomograms(self.extractedParts, relionParticles,
                                          expectedSetSize=expectedSetSize,
                                          expectedSRate=expectedSRate,
                                          expectedBoxSize=expectedBoxSize,
                                          orientedParticles=True,
                                          isStack2d=are2dStacks,
                                          noImgs=self.nImgs,
                                          expectedExtension='.mrcs')  # 2d stacks

        # Check the classes
        self.checkClasses(outClasses,
                          expectedSetSize=self.nClasses,
                          expectedSRate=expectedSRate)


class TestRelion5TomoRefine(TestRelion5RefineCycleBase):

    @classmethod
    def _runPreviousProtocols(cls):
        super()._runPreviousProtocols()
        cls.importedRef = cls._runImportReference(filePath=cls.ds.getFile(DataSetRe4STATuto.recParticleBin6.name),
                                                  samplingRate=cls.unbinnedSRate * cls.binFactor6)
        cls.protExtract = cls._runExtractSubtomos(inParticles=cls.importedCoords,
                                                  inputCtfTs=cls.importedCtfs,
                                                  inputTS=cls.tsWithAlignment,
                                                  binningFactor=cls.binFactor6,
                                                  boxSize=cls.boxSizeBin6,
                                                  croppedBoxSize=cls.croppedBoxSizeBin6,
                                                  gen2dParticles=True)
        cls.extractedParts = getattr(cls.protExtract, cls.protExtract._possibleOutputs.relionParticles.name, None)

    def testAutoRefine(self):
        if IS_RE_50:
            binningFactor = self.binFactor6
            protAutoRefine = self._runAutoRefine()
            are2dStacks = True
            relionParticles = getattr(protAutoRefine, protAutoRefine._possibleOutputs.relionParticles.name, None)
            # Check RelionTomoMetadata: only the particles file is generated
            self._checkRe4Metadata(relionParticles,
                                   noParticles=self.nParticles,
                                   tomogramsFile=self.protExtract._getExtraPath(IN_TOMOS_STAR),
                                   particlesFile=protAutoRefine._getExtraPath('_data.star'),
                                   trajectoriesFile=None,
                                   referenceFscFile=None,
                                   relionBinning=binningFactor,
                                   are2dStacks=are2dStacks,
                                   tsSamplingRate=self.unbinnedSRate)

            # Check the set of pseudosubtomograms
            expectedSetSize = self.nParticles
            expectedSRate = self.unbinnedSRate * binningFactor
            expectedBoxSize = self.croppedBoxSizeBin6

            self.checkRefinedSubtomograms(self.extractedParts, relionParticles,
                                          expectedSetSize=expectedSetSize,
                                          expectedSRate=expectedSRate,
                                          expectedBoxSize=expectedBoxSize,
                                          orientedParticles=True,
                                          isStack2d=are2dStacks,
                                          noImgs=self.nImgs,
                                          expectedExtension='.mrcs')  # 2d stacks

            # Check the output volume
            avg = getattr(protAutoRefine, protAutoRefine._possibleOutputs.average.name, None)
            self.checkAverage(avg,
                              expectedSRate=expectedSRate,
                              expectedBoxSize=expectedBoxSize,
                              hasHalves=True)
        else:
            print(yellowStr('Relion 4 detected. Test for protocol "Auto-refine" skipped.'))

    @classmethod
    def _runAutoRefine(cls):
        print(magentaStr("\n==> Refining the particles:"))
        protAutoRefine = cls.newProtocol(ProtRelionRefineSubtomograms,
                                         inReParticles=cls.extractedParts,
                                         referenceVolume=cls.importedRef,
                                         isMapAbsoluteGreyScale=False,
                                         initialLowPassFilterA=60,
                                         symmetry=cls.symmetry,
                                         maskDiameter=500,
                                         useFinerAngularSampling=False,
                                         priorWidthTiltAngle=10,
                                         pooledSubtomos=6,
                                         doGpu=True,
                                         gpusToUse='0',
                                         numberOfMpi=3,
                                         binThreads=3)
        cls.launchProtocol(protAutoRefine)
        return protAutoRefine


class TestRelion5PostProcess(TestRelion5RefineCycleBase):
    boxSizeBin4 = DataSetRe4STATuto.boxSizeBin4.value
    croppedBoxSizeBin4 = DataSetRe4STATuto.croppedBoxSizeBin4.value

    @classmethod
    def _runPreviousProtocols(cls):
        cls.importedFscMask = cls._runImportReference(filePath=cls.ds.getFile(DataSetRe4STATuto.maskFscBin2.name),
                                                      samplingRate=cls.unbinnedSRate * cls.binFactor2)

    def test_runPostProcess(self):
        if IS_RE_50:
            sRateBin2 = self.unbinnedSRate * self.binFactor2
            print(magentaStr("\n==> Importing the reference volume with halves:"))
            protImportRef = self.newProtocol(ProtImportVolumes,
                                             filesPath=self.ds.getFile(DataSetRe4STATuto.recParticleBin2.name),
                                             samplingRate=sRateBin2,
                                             setHalfMaps=True,
                                             half1map=self.ds.getFile(DataSetRe4STATuto.recParticleHalf1Bin2.name),
                                             half2map=self.ds.getFile(DataSetRe4STATuto.recParticleHalf2Bin2.name))
            self.launchProtocol(protImportRef)
            refWithHalves = getattr(protImportRef, ProtImportVolumes._possibleOutputs.outputVolume.name, None)
            print(magentaStr("\n==> Post-processing:"))
            protPostProcess = self.newProtocol(ProtRelionPostProcess,
                                               inVolume=refWithHalves,
                                               solventMask=self.importedFscMask)
            self.launchProtocol(protPostProcess)

            # Check the output volume
            avg = getattr(protPostProcess, protPostProcess._possibleOutputs.postProcessVolume.name, None)
            self.checkAverage(avg,
                              expectedSRate=sRateBin2,
                              expectedBoxSize=self.croppedBoxSizeBin2,
                              hasHalves=False)

            # Check the FSC (it must contain an extended attribute pointing to the postprocess star file that can be
            # introduced as input in the CTF refinement and bayesian polishing
            fsc = getattr(protPostProcess, protPostProcess._possibleOutputs.outputFSC.name, None)
            self.assertTrue(hasattr(fsc, POSTPROCESS_STAR_FIELD))
            self.assertTrue(exists(getattr(fsc, POSTPROCESS_STAR_FIELD).get()))
        else:
            print(yellowStr('Relion 4 detected. Test for protocol "Post-process" skipped.'))
