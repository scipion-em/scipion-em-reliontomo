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
from pwem.emlib.image import ImageHandler
from pwem.protocols import ProtImportMask, ProtImportVolumes
from pwem.protocols.protocol_import.masks import ImportMaskOutput
from pyworkflow.tests import setupTestProject, DataSet
from pyworkflow.utils import magentaStr
from reliontomo import Plugin
from reliontomo.constants import OUT_PARTICLES_STAR, IN_TOMOS_STAR, FILE_NOT_FOUND
from reliontomo.convert.convertBase import getTransformInfoFromCoordOrSubtomo
from reliontomo.protocols import ProtImportCoordinates3DFromStar, ProtRelion5ExtractSubtomos, \
    ProtRelion5ReconstructParticle, ProtRelionDeNovoInitialModel
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


class TestRelion5RefineCycleBase(TestBaseCentralizedLayer):
    ds = None
    tsWithAlignment = None
    importedTs = None
    importedCtfs = None
    importedTomos = None
    importedCoords = None
    extractedParticles = None
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
    def _runExtractSubtomos(cls, inParticles=None, binningFactor=None, boxSize=None, croppedBoxSize=None,
                            gen2dParticles=False):
        extractMsg = 'Extract' if type(inParticles) is SetOfCoordinates3D else 'Re-extract'
        partTypeMsg = '2D' if gen2dParticles else '3D'
        print(magentaStr(f"\n==> {extractMsg}ing the particles from the tilt-series:\n"
                         f"\t- Particles type selected: {partTypeMsg}"))
        protExtract = cls.newProtocol(ProtRelion5ExtractSubtomos,
                                      inputCtfTs=cls.importedCtfs,
                                      inReParticles=inParticles,
                                      inputTS=cls.tsWithAlignment,
                                      binningFactor=binningFactor,
                                      boxSize=boxSize,
                                      croppedBoxSize=croppedBoxSize,
                                      write2dStacks=gen2dParticles,
                                      numberOfMpi=2,  # There are 2 tomograms in the current dataset
                                      numberOfThreads=4)
        protExtract.setObjLabel(f"{extractMsg} subtomos {partTypeMsg}")
        return cls.launchProtocol(protExtract)

    @classmethod
    def _runImportReference(cls):
        print(magentaStr("\n==> Importing the reference volume:"))
        protImportRef = cls.newProtocol(ProtImportVolumes,
                                        filesPath=cls.ds.getFile(DataSetRe4STATuto.initModelRelion.name),
                                        samplingRate=DataSetRe4STATuto.sRateBin4.value)
        cls.launchProtocol(protImportRef)
        return getattr(protImportRef, ProtImportVolumes._possibleOutputs.outputVolume.name, None)

    def _checkRe4Metadata(self, pSubtomoSet, tomogramsFile=None, particlesFile=None, trajectoriesFile=None,
                          manifoldsFile=None, referenceFscFile=None, relionBinning=None):
        self.assertEqual(DataSetRe4STATuto.nCoordsTotal.value, pSubtomoSet.getNReParticles())
        self.assertEqual(self.unbinnedSRate, pSubtomoSet.getTsSamplingRate())
        self.assertEqual(tomogramsFile, pSubtomoSet.getTomogramsStar())
        self.assertEqual(particlesFile, pSubtomoSet.getParticlesStar())
        self.assertEqual(trajectoriesFile, pSubtomoSet.getTrajectoriesStar())
        self.assertEqual(manifoldsFile, pSubtomoSet.getManifolds())
        self.assertEqual(referenceFscFile, pSubtomoSet.getReferenceFsc())
        self.assertEqual(relionBinning, pSubtomoSet.getRelionBinning())

    def _checkPseudosubtomograms(self, inCoords, outSubtomos, expectedSetSize=-1, expectedSRate=-1, expectedBoxSize=-1,
                                 convention=TR_SCIPION, orientedParticles=False, are2dStacks=None, nTiltImages=None,
                                 expectedRelionBinning=None, tsSamplingRate=None, expectedNParticles=None):
        # This way the box size check will be skipped in checkExtractedSubtomos and the corresponding box size will
        # be checked below here
        expBoxSize = None if are2dStacks else expectedBoxSize
        # Check the SubTomogram part, as Relionparticles inherit from them
        self.checkExtractedSubtomos(inCoords, outSubtomos,
                                    expectedSetSize=expectedSetSize,
                                    expectedSRate=expectedSRate,
                                    expectedBoxSize=expBoxSize,
                                    convention=convention,
                                    orientedParticles=orientedParticles)

        # Check some remaining specific relionParticle attributes
        self.assertEqual(outSubtomos.getRelionBinning(), expectedRelionBinning)
        self.assertEqual(outSubtomos.getTsSamplingRate(), tsSamplingRate)
        self.assertEqual(outSubtomos.getCurrentSamplingRate(), tsSamplingRate * expectedRelionBinning)
        self.assertEqual(outSubtomos.getNReParticles(), expectedNParticles)
        self.assertEqual(outSubtomos.areRe5Particles(), True if Plugin.isRe50() else False)

        if are2dStacks:
            self.assertTrue(outSubtomos.are2dStacks())
            ih = ImageHandler()
            for pSubtomo in outSubtomos:
                x, y, z, _ = ih.getDimensions(pSubtomo.getFileName())
                self.assertEqual(x, expectedBoxSize)
                self.assertEqual(y, expectedBoxSize)
                self.assertLessEqual(z, nTiltImages)
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

    def _runTestExtractSubtomos(self, inParticles=None, binningFactor=None, boxSize=None, croppedBoxSize=None,
                                are2dParticles=False, nTiltImages=None):
        protExtract = self._runExtractSubtomos(inParticles=inParticles,
                                               binningFactor=binningFactor,
                                               boxSize=boxSize,
                                               croppedBoxSize=croppedBoxSize,
                                               gen2dParticles=are2dParticles)
        outParticles = getattr(protExtract, ProtRelion5ExtractSubtomos._possibleOutputs.relionParticles.name, None)
        # Check RelionTomoMetadata: both particles and tomograms files are generated
        self._checkRe4Metadata(outParticles,
                               tomogramsFile=protExtract._getExtraPath(IN_TOMOS_STAR),
                               particlesFile=protExtract._getExtraPath(OUT_PARTICLES_STAR),
                               trajectoriesFile=None,
                               manifoldsFile=None,
                               referenceFscFile=None,
                               relionBinning=binningFactor)
        # Check that the projected coordinates have been generated
        self.assertIsNotNone(
            getattr(protExtract, ProtRelion5ExtractSubtomos._possibleOutputs.projected2DCoordinates.name, None),
            msg='The projected coordinates were not generated.')
        # Check the pseudo-subtomograms
        unbinnedPixSize = self.unbinnedSRate
        nParticles = DataSetRe4STATuto.nCoordsTotal.value
        self._checkPseudosubtomograms(self.importedCoords, outParticles,
                                      expectedSetSize=nParticles,
                                      expectedSRate=unbinnedPixSize * binningFactor,
                                      expectedBoxSize=croppedBoxSize,
                                      orientedParticles=True,  # The imported coordinates are oriented
                                      are2dStacks=are2dParticles,
                                      nTiltImages=nTiltImages,
                                      expectedRelionBinning=binningFactor,
                                      tsSamplingRate=unbinnedPixSize,
                                      expectedNParticles=nParticles)
        return outParticles

    def testExtractSubtomos3d(self):
        self._runTestExtractSubtomos(inParticles=self.importedCoords,
                                     binningFactor=self.binFactor6,
                                     boxSize=self.boxSizeBin6,
                                     croppedBoxSize=self.croppedBoxSizeBin6,
                                     are2dParticles=False)

    def testExtractSubtomos2d(self):
        extractedParticles = self._runTestExtractSubtomos(inParticles=self.importedCoords,
                                                          binningFactor=self.binFactor6,
                                                          boxSize=self.boxSizeBin6,
                                                          croppedBoxSize=self.croppedBoxSizeBin6,
                                                          are2dParticles=True,
                                                          nTiltImages=40)
        # Run a re-extraction
        self._runTestExtractSubtomos(inParticles=extractedParticles,
                                     binningFactor=self.binFactor2,
                                     boxSize=self.boxSizeBin2,
                                     croppedBoxSize=self.croppedBoxSizeBin2,
                                     are2dParticles=True,
                                     nTiltImages=40)


class TestRelionTomoRecParticleFromTs(TestRelion5RefineCycleBase):

    def testRecParticleFromTS_01(self):
        self._runTest(are2dStacks=True, hasHalves=False)  # Halves are generated after having run a refinement

    def testRecParticleFromTS_02(self):
        self._runTest(are2dStacks=False, hasHalves=False)

    def _runTest(self, are2dStacks=True, hasHalves=None):
        binningFactor = self.binFactor6
        bozSize = self.boxSizeBin6
        croppedBoxSize = self.croppedBoxSizeBin6
        protExtract = self._runExtractSubtomos(inParticles=self.importedCoords,
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
                                            numberOfThreads=4)
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


class TestRelionTomoEditStar(TestRelion5RefineCycleBase):
    editStarOperationDict = genEnumParamDict(OPERATION_LABELS)
    editStarLabelsDict = genEnumParamDict(LABELS_TO_OPERATE_WITH)
    editTestsTol = 0.01

    @classmethod
    def _runPreviousProtocols(cls):
        super()._runPreviousProtocols()
        cls.importedRef = cls._runImportReference()
        cls.protExtract = cls._runExtractSubtomos(inParticles=cls.importedCoords,
                                                  binningFactor=cls.binFactor6,
                                                  boxSize=cls.boxSizeBin6,
                                                  croppedBoxSize=cls.croppedBoxSizeBin6,
                                                  gen2dParticles=True)
        cls.inReParticles = getattr(cls.protExtract, cls.protExtract._possibleOutputs.relionParticles.name, None)

    def testEditStar_shiftCenter(self):
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

    def testEditStar_removeDuplicates(self):
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
        self._checkPseudosubtomograms(self.importedCoords, outParticles,
                                      expectedSetSize=expectedSetSize,
                                      expectedSRate=unbinnedPixSize * binningFactor,
                                      expectedBoxSize=self.croppedBoxSizeBin6,
                                      orientedParticles=True,  # The imported coordinates are oriented
                                      are2dStacks=True,
                                      nTiltImages=outParticles.are2dStacks(),
                                      expectedRelionBinning=binningFactor,
                                      tsSamplingRate=unbinnedPixSize,
                                      expectedNParticles=expectedSetSize)

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


class TestRelionTomoGenInitialModel(TestRelion5RefineCycleBase):

    def testInitialModel(self):
        recVol = self._genInitialModel()
        self.checkAverage(recVol,
                          expectedSRate=self.unbinnedSRate * self.binFactor6,
                          expectedBoxSize=self.croppedBoxSizeBin6,
                          hasHalves=False)

    @classmethod
    def _genInitialModel(cls, are2dStacks=True):
        protExtract = cls._runExtractSubtomos(
            inParticles=cls.importedCoords,
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
                                           numberOfThreads=6)
        cls.launchProtocol(protInitialModel)
        return getattr(protInitialModel, ProtRelionDeNovoInitialModel._possibleOutputs.average.name, None)


class TestRelionTomo3dClassify(TestRelion5RefineCycleBase):
    nClasses = 2

    @classmethod
    def _runPreviousProtocols(cls):
        super()._runPreviousProtocols()
        cls.importedRef = cls._runImportReference()
        protExtract = cls._runExtractSubtomos(inParticles=cls.importedCoords,
                                              binningFactor=cls.binFactor6,
                                              boxSize=cls.boxSizeBin6,
                                              croppedBoxSize=cls.croppedBoxSizeBin6,
                                              gen2dParticles=True)
        cls.extractedPars = getattr(protExtract, ProtRelion5ExtractSubtomos._possibleOutputs.relionParticles.name, None)
#
#     def testCl3d(self):
#         protCl3d = self._run3dClassify()
#         self._checkCl3dResults(protCl3d)
#
#     def testCl3dWithAlignment(self):
#         protCl3d = self._run3dClassify(doAlingment=True)
#         self._checkCl3dResults(protCl3d, onlyClassify=False)
#
    # @classmethod
    # def _run3dClassify(cls, doAlingment=False):
    #     paramsDict = {'inReParticles': cls.inReParticles,
    #                   'referenceVolume': cls.importedRef,
    #                   'numberOfClasses': cls.nClasses,
    #                   'initialLowPassFilterA': 30,
    #                   'symmetry': cls.symmetry,
    #                   'maskDiameter': 230,
    #                   'pooledSubtomos': 6,
    #                   'numberOfMpi': 1,
    #                   'numberOfThreads': 3}
    #
    #     if doAlingment:
    #         paramsDict['doImageAlignment'] = True
    #         paramsDict['nIterations'] = 1  # Prevent error "No orientation was found as better than any other"
    #         # because of the test dataset
    #         paramsDict['doGpu'] = True
    #         paramsDict['gpusToUse'] = '0'
    #         label = 'cl3d with alignment'
    #         print(magentaStr("\n==> Classifying and aligning the psudosubtomograms:"))
    #     else:
    #         paramsDict['nIterations'] = 3
    #         label = 'cl3d'
    #         print(magentaStr("\n==> Classifying the psudosubtomograms:"))
    #
    #     protCl3d = cls.newProtocol(ProtRelion3DClassifySubtomograms, **paramsDict)
    #     protCl3d.setObjLabel(label)
    #     cls.launchProtocol(protCl3d)
    #     return protCl3d
#
#     def _checkCl3dResults(self, protCl3d, onlyClassify=True):
#         relionParticles = getattr(protCl3d, ProtRelion3DClassifySubtomograms._possibleOutputs.relionParticles.name,
#                                   None)
#         outClasses = getattr(protCl3d, ProtRelion3DClassifySubtomograms._possibleOutputs.classes.name, None)
#         # Check RelionTomoMetadata: only the particles file is generated
#         self._checkRe4Metadata(relionParticles,
#                                tomogramsFile=self.protExtract._getExtraPath(OUT_TOMOS_STAR),
#                                particlesFile=protCl3d._getExtraPath(OUT_PARTICLES_STAR),
#                                trajectoriesFile=None,
#                                manifoldsFile=None,
#                                referenceFscFile=None,
#                                relionBinning=self.binFactor4
#                                )
#         # Check the set of pseudosubtomograms
#         expectedSetSize = DataSetRe4STATuto.nCoordsTotal.value
#         expectedSRate = self.unbinnedSRate * self.binFactor4
#         expectedBoxSize = self.croppedBoxSizeBin6
#         if onlyClassify:
#             inCoords = self.protExtract.inputCoords.get()
#             self._checkPseudosubtomograms(inCoords, relionParticles,
#                                           expectedSetSize=expectedSetSize,
#                                           expectedSRate=expectedSRate,
#                                           expectedBoxSize=expectedBoxSize,
#                                           orientedParticles=True)  # They wer picked with PySeg
#         else:
#             self.checkRefinedSubtomograms(self.inReParticles, relionParticles,
#                                           expectedSetSize=expectedSetSize,
#                                           expectedSRate=expectedSRate,
#                                           expectedBoxSize=expectedBoxSize,
#                                           orientedParticles=True)  # They wer picked with PySeg
#
#         # Check the classes
#         self.checkClasses(outClasses,
#                           expectedSetSize=self.nClasses,
#                           expectedSRate=expectedSRate)


# class TestRelionTomoRefine(TestRelion5RefineCycleBase):
#
#     @classmethod
#     def _runPreviousProtocols(cls):
#         super()._runPreviousProtocols()
#         cls.importedRef = cls._runImportReference()
#         cls.protExtract = cls._runExtractSubtomos()
#         protMakePSubtomos = cls._makePSubtomograms()
#         cls.inReParticles = getattr(protMakePSubtomos,
#                                     ProtRelionMakePseudoSubtomograms._possibleOutputs.relionParticles.name, None)
#
#     def testAutoRefine(self):
#         protAutoRefine = self._runAutoRefine()
#         relionParticles = getattr(protAutoRefine, ProtRelionRefineSubtomograms._possibleOutputs.relionParticles.name,
#                                   None)
#         # Check RelionTomoMetadata: only the particles file is generated
#         self._checkRe4Metadata(relionParticles,
#                                tomogramsFile=self.protExtract._getExtraPath(OUT_TOMOS_STAR),
#                                particlesFile=protAutoRefine._getExtraPath('_data.star'),
#                                trajectoriesFile=None,
#                                manifoldsFile=None,
#                                referenceFscFile=None,
#                                relionBinning=self.binFactor4)
#
#         # Check the set of pseudosubtomograms
#         expectedSetSize = DataSetRe4STATuto.nCoordsTotal.value
#         expectedSRate = self.unbinnedSRate * self.binFactor4
#         expectedBoxSize = self.croppedBoxSizeBin6
#
#         self.checkRefinedSubtomograms(self.inReParticles, relionParticles,
#                                       expectedSetSize=expectedSetSize,
#                                       expectedSRate=expectedSRate,
#                                       expectedBoxSize=expectedBoxSize,
#                                       orientedParticles=True)
#
#         # Check the output volume
#         avg = getattr(protAutoRefine, ProtRelionRefineSubtomograms._possibleOutputs.average.name, None)
#         self.checkAverage(avg,
#                           expectedSRate=expectedSRate,
#                           expectedBoxSize=expectedBoxSize,
#                           hasHalves=True)
#
#     @classmethod
#     def _runAutoRefine(cls):
#         print(magentaStr("\n==> Refining the particles:"))
#         protAutoRefine = cls.newProtocol(ProtRelionRefineSubtomograms,
#                                          inReParticles=cls.inReParticles,
#                                          referenceVolume=cls.importedRef,
#                                          isMapAbsoluteGreyScale=False,
#                                          initialLowPassFilterA=50,
#                                          symmetry=cls.symmetry,
#                                          maskDiameter=230,
#                                          useFinerAngularSampling=True,
#                                          pooledSubtomos=6,
#                                          doGpu=True,
#                                          gpusToUse='0',
#                                          numberOfMpi=3,
#                                          numberOfThreads=3)
#         cls.launchProtocol(protAutoRefine)
#         return protAutoRefine
