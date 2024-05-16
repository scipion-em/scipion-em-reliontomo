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
from pwem.convert.transformations import translation_from_matrix
from pwem.emlib.image import ImageHandler
from pwem.protocols import ProtImportMask, ProtImportVolumes
from pwem.protocols.protocol_import.masks import ImportMaskOutput
from pyworkflow.tests import setupTestProject, DataSet
from pyworkflow.utils import magentaStr
from reliontomo.constants import OUT_TOMOS_STAR, OUT_PARTICLES_STAR, IN_PARTICLES_STAR, POSTPROCESS_DIR, \
    POST_PROCESS_MRC, IN_TOMOS_STAR, FILE_NOT_FOUND
from reliontomo.convert.convertBase import getTransformInfoFromCoordOrSubtomo
from reliontomo.protocols import ProtImportCoordinates3DFromStar, ProtRelion5ExtractSubtomos, \
    ProtRelionMakePseudoSubtomograms, ProtRelionDeNovoInitialModel, ProtRelionRefineSubtomograms, \
    ProtRelionReconstructParticle, ProtRelionTomoReconstruct, ProtRelionEditParticlesStar
from reliontomo.protocols.protocol_3d_classify_subtomograms import ProtRelion3DClassifySubtomograms
from reliontomo.protocols.protocol_edit_particles_star import OPERATION_LABELS, LABELS_TO_OPERATE_WITH, ANGLES, \
    OP_ADDITION, OP_MULTIPLICATION, COORDINATES, OP_SET_TO
from reliontomo.protocols.protocol_edit_particles_star import outputObjects as editStarOutputs
from reliontomo.protocols.protocol_rec_tomogram import SINGLE_TOMO, ALL_TOMOS
from reliontomo.tests import DataSetRe4Tomo
from reliontomo.utils import genEnumParamDict
from tomo.constants import TR_SCIPION
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
    croppedBoxSizeBin6 = DataSetRe4STATuto.croppedBoxSizeBin4.value
    boxSizeBin6 = DataSetRe4STATuto.boxSizeBin4.value
    tsIds = ['TS_03', 'TS_54']
    binFactor4 = 4
    binFactor6 = 6
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
                                       samplingRate=DataSetRe4STATuto.unbinnedPixSize.value,
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
                                         samplingRate=DataSetRe4STATuto.unbinnedPixSize.value)
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
                                           samplingRate=DataSetRe4STATuto.unbinnedPixSize.value)
        cls.launchProtocol(protImportCoords)
        outCoords = getattr(protImportCoords, ProtImportCoordinates3DFromStar._possibleOutputs.coordinates.name, None)
        return outCoords

    @classmethod
    def _runExtractSubtomos(cls, binningFactor, boxSize, croppedBoxSize, gen2dParticles=False):
        partTypeMsg = '2D' if gen2dParticles else '3D'
        print(magentaStr(f"\n==> Extracting the particles from the tilt-series:\n"
                         f"\t- Particles type selected: {partTypeMsg}"))
        protExtract = cls.newProtocol(ProtRelion5ExtractSubtomos,
                                      inputCtfTs=cls.importedCtfs,
                                      inReParticles=cls.importedCoords,
                                      inputTS=cls.tsWithAlignment,
                                      binningFactor=binningFactor,
                                      boxSize=boxSize,
                                      croppedBoxSize=croppedBoxSize,
                                      write2dStacks=gen2dParticles,
                                      numberOfMpi=2,  # There are 2 tomograms in the current dataset
                                      numberOfThreads=3)
        protExtract.setObjLabel(f"Extract subtomos {partTypeMsg}")
        return cls.launchProtocol(protExtract)

    # @classmethodb
    # def _makePSubtomograms(cls):
    #     print(magentaStr("\n==> Making the psudosubtomograms:"))
    #     inRelionParticles = getattr(cls.protExtract, ProtRelion5ExtractSubtomos._possibleOutputs.relionParticles.name, None)
    #     protMakePsubtomos = cls.newProtocol(ProtRelionMakePseudoSubtomograms,
    #                                         inReParticles=inRelionParticles,
    #                                         boxSize=cls.boxSizeBin6,
    #                                         croppedBoxSize=cls.croppedBoxSizeBin6,
    #                                         binningFactor=cls.binFactor4,
    #                                         outputInFloat16=False,
    #                                         numberOfThreads=5,
    #                                         numberOfMpi=3)
    #     return cls.launchProtocol(protMakePsubtomos)

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
        unbinnedSRate = DataSetRe4STATuto.unbinnedPixSize.value
        self.assertEqual(DataSetRe4STATuto.nCoordsTotal.value, pSubtomoSet.getNReParticles())
        self.assertEqual(unbinnedSRate, pSubtomoSet.getTsSamplingRate())
        self.assertEqual(tomogramsFile, pSubtomoSet.getTomogramsStar())
        self.assertEqual(particlesFile, pSubtomoSet.getParticlesStar())
        self.assertEqual(trajectoriesFile, pSubtomoSet.getTrajectoriesStar())
        self.assertEqual(manifoldsFile, pSubtomoSet.getManifolds())
        self.assertEqual(referenceFscFile, pSubtomoSet.getReferenceFsc())
        self.assertEqual(relionBinning, pSubtomoSet.getRelionBinning())

    def _checkPseudosubtomograms(self, inCoords, outSubtomos, expectedSetSize=-1, expectedSRate=-1, expectedBoxSize=-1,
                                 convention=TR_SCIPION, orientedParticles=False, are2dStacks=None, nTiltImages=None):
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
        if are2dStacks:
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
            for pSubtomo in outSubtomos:
                self.assertEqual(len(pSubtomo.getVisibleFrames()), 1)
                self.assertTrue(exists(pSubtomo.getCtfFile()))
                self.assertTrue(pSubtomo.getTsId() in self.tsIds)


class TestRelion5TomoExtractSubtomos(TestRelion5RefineCycleBase):

    def _runTestExtractSubtomos(self, are2dParticles=False, nTiltImages=None):
        protExtract = self._runExtractSubtomos(self.binFactor6, self.boxSizeBin6, self.croppedBoxSizeBin6,
                                               gen2dParticles=are2dParticles)
        outParticles = getattr(protExtract, ProtRelion5ExtractSubtomos._possibleOutputs.relionParticles.name, None)
        # Check RelionTomoMetadata: both particles and tomograms files are generated
        self._checkRe4Metadata(outParticles,
                               tomogramsFile=protExtract._getExtraPath(IN_TOMOS_STAR),
                               particlesFile=protExtract._getExtraPath(OUT_PARTICLES_STAR),
                               trajectoriesFile=None,
                               manifoldsFile=None,
                               referenceFscFile=None,
                               relionBinning=self.binFactor6)
        # Check that the projected coordinates have been generated
        self.assertIsNotNone(
            getattr(protExtract, ProtRelion5ExtractSubtomos._possibleOutputs.projected2DCoordinates.name, None),
            msg='The projected coordinates were not generated.')
        # Check the pseudo-subtomograms
        self._checkPseudosubtomograms(self.importedCoords, outParticles,
                                      expectedSetSize=DataSetRe4STATuto.nCoordsTotal.value,
                                      expectedSRate=DataSetRe4STATuto.unbinnedPixSize.value * self.binFactor6,
                                      expectedBoxSize=self.croppedBoxSizeBin6,
                                      orientedParticles=True,  # The imported coordinates are oriented
                                      are2dStacks=are2dParticles,
                                      nTiltImages=nTiltImages)

    def testExtractSubtomos3d(self):
        self._runTestExtractSubtomos()

    def testExtractSubtomos2d(self):
        self._runTestExtractSubtomos(are2dParticles=True, nTiltImages=40)


class TestRelionTomoRecTomograms(TestRelion5RefineCycleBase):
    binFactor = 8
    expectedTomoDims = [464, 464, 140]

    @classmethod
    def _runPreviousProtocols(cls):
        super()._runPreviousProtocols()
        cls.protExtract = cls._runExtractSubtomos()

    def testRecSingleTomoFromPrep(self):
        expectedSize = 1
        tomoSet = self._runRecFromPrepare_SingleTomo()
        self._checkTomograms(tomoSet, expectedSetSize=expectedSize)

    def testRecAllTomosFromPrep(self):
        expectedSize = 2
        tomoSet = self._runRecFromPrepare_AllTomos()
        self._checkTomograms(tomoSet, expectedSetSize=expectedSize)

    @classmethod
    def _runRecFromPrepare_SingleTomo(cls):
        print(magentaStr("\n==> Reconstructing one of the tomograms with Relion:"))
        protRelionRec = cls._recFromPrepare(SINGLE_TOMO, tomoId='TS_54')  #
        return protRelionRec

    @classmethod
    def _runRecFromPrepare_AllTomos(cls):
        print(magentaStr("\n==> Reconstructing all the tomograms with Relion:"))
        protRelionRec = cls._recFromPrepare(ALL_TOMOS)
        return protRelionRec

    @classmethod
    def _recFromPrepare(cls, recMode, tomoId=None):
        protRelionRec = cls.newProtocol(ProtRelionTomoReconstruct,
                                        protExtract=cls.protExtract,
                                        recTomoMode=recMode,
                                        tomoId=tomoId,
                                        binFactor=cls.binFactor)
        cls.launchProtocol(protRelionRec)
        return getattr(protRelionRec, ProtRelionTomoReconstruct._possibleOutputs.tomograms.name, None)

    def _checkTomograms(self, inTomoSet, expectedSetSize=None):
        self.checkTomograms(inTomoSet,
                            expectedSetSize=expectedSetSize,
                            expectedSRate=DataSetRe4STATuto.unbinnedPixSize.value * self.binFactor,
                            expectedDimensions=self.expectedTomoDims)


class TestRelionTomoRecParticleFromTs(TestRelion5RefineCycleBase):
    binFactor2 = 2
    croppedBoxSizeBin2 = DataSetRe4STATuto.croppedBoxSizeBin2.value
    boxSizeBin2 = DataSetRe4STATuto.boxSizeBin2.value

    @classmethod
    def _runPreviousProtocols(cls):
        super()._runPreviousProtocols()
        cls.protExtract = cls._runExtractSubtomos()

    def testRecParticleFromTS(self):
        protRecPartFromTS = self._recParticleFromTS()
        # Check the results
        self._checkResults(protRecPartFromTS)

    def testRecParticleFromTsWithSolventMask(self):
        protRecPartFromTsWithSolvent = self._recParticleFromTsWithSolventMask()
        # Check the results
        fscFile = protRecPartFromTsWithSolvent._getExtraPath(POSTPROCESS_DIR, POST_PROCESS_MRC.replace('.mrc', '.star'))
        self._checkResults(protRecPartFromTsWithSolvent, fscFile=fscFile)

    @classmethod
    def _genRecPartFromTsDict(cls, solventMask=None):
        inReParticles = getattr(cls.protExtract, ProtRelion5ExtractSubtomos._possibleOutputs.relionParticles.name, None)
        paramsDict = {'inReParticles': inReParticles,
                      'boxSize': cls.boxSizeBin2,
                      'croppedBoxSize': cls.croppedBoxSizeBin2,
                      'binningFactor': cls.binFactor2,
                      'symmetry': cls.symmetry,
                      'outputInFloat16': False}
        if solventMask:
            paramsDict['solventMask'] = solventMask
        return paramsDict

    @classmethod
    def _recParticleFromTS(cls):
        print(magentaStr("\n==> Reconstructing the particle from the TS using a binning factor of 2:"))
        protRecPartFromTS = cls.newProtocol(ProtRelionReconstructParticle, **cls._genRecPartFromTsDict())
        protRecPartFromTS.setObjLabel('rec part from ts')
        cls.launchProtocol(protRecPartFromTS)
        return protRecPartFromTS

    @classmethod
    def _recParticleFromTsWithSolventMask(cls):
        print(magentaStr("\n==> Reconstructing the particle from the TS using a binning factor of 2 "
                         "with solvent mask:"))
        solventMaskBin2 = cls._runImportFscMask()
        paramsDict = cls._genRecPartFromTsDict(solventMask=solventMaskBin2)
        protRecPartFromTS = cls.newProtocol(ProtRelionReconstructParticle, **paramsDict)
        protRecPartFromTS.setObjLabel('rec part from ts with solvent mask')
        cls.launchProtocol(protRecPartFromTS)
        return protRecPartFromTS

    def _checkResults(self, protRecPartFromTS, fscFile=None):
        reParticles = getattr(protRecPartFromTS, ProtRelionReconstructParticle._possibleOutputs.relionParticles.name,
                              None)
        recVol = getattr(protRecPartFromTS, ProtRelionReconstructParticle._possibleOutputs.average.name, None)
        # Check the metadata
        self._checkRe4Metadata(reParticles,
                               tomogramsFile=self.protExtract._getExtraPath(OUT_TOMOS_STAR),
                               particlesFile=protRecPartFromTS._getExtraPath(IN_PARTICLES_STAR),
                               trajectoriesFile=None,
                               manifoldsFile=None,
                               referenceFscFile=fscFile,
                               relionBinning=self.binFactor2)
        # Check the reconstructed volume
        self.checkAverage(recVol,
                          expectedSRate=DataSetRe4STATuto.unbinnedPixSize.value * self.binFactor2,
                          expectedBoxSize=self.croppedBoxSizeBin2,
                          hasHalves=True)


# class TestRelionTomoMakePseudoSubtomos(TestRelion5RefineCycleBase):
#
#     @classmethod
#     def _runPreviousProtocols(cls):
#         super()._runPreviousProtocols()
#         cls.protExtract = cls._runExtractSubtomos()
#
#     def testMakePSubtomos(self):
#         protMakePSubtomos = self._makePSubtomograms()
#         relionParticles = getattr(protMakePSubtomos,
#                                   ProtRelionMakePseudoSubtomograms._possibleOutputs.relionParticles.name, None)
#         # Check RelionTomoMetadata: only the particles file is generated
#         self._checkRe4Metadata(relionParticles,
#                                tomogramsFile=self.protExtract._getExtraPath(OUT_TOMOS_STAR),
#                                particlesFile=protMakePSubtomos._getExtraPath(OUT_PARTICLES_STAR),
#                                trajectoriesFile=None,
#                                manifoldsFile=None,
#                                referenceFscFile=None,
#                                relionBinning=self.binFactor4
#                                )
#         # Check the set of pseudosubtomograms
#         inCoords = self.protExtract.inputCoords.get()
#         self._checkPseudosubtomograms(inCoords, relionParticles,
#                                       expectedSetSize=DataSetRe4STATuto.nCoordsTotal.value,
#                                       expectedSRate=DataSetRe4STATuto.unbinnedPixSize.value * self.binFactor4,
#                                       expectedBoxSize=self.croppedBoxSizeBin6,
#                                       orientedParticles=True)  # The imported coordinates were picked using pyseg, so
#         # they are oriented


class TestRelionTomoEditStar(TestRelion5RefineCycleBase):
    editStarOperationDict = genEnumParamDict(OPERATION_LABELS)
    editStarLabelsDict = genEnumParamDict(LABELS_TO_OPERATE_WITH)
    editTestsTol = 0.01

    @classmethod
    def _runPreviousProtocols(cls):
        super()._runPreviousProtocols()
        cls.importedRef = cls._runImportReference()
        cls.protExtract = cls._runExtractSubtomos()
        cls.inReParticles = getattr(cls.protExtract, cls.protExtract._possibleOutputs.relionParticles.name, None)

    def testEditStar_shiftCenter(self):
        # Values edited: shiftX = 4, shiftY = 2, shiftZ = 3
        protEdit = self._editStar_shiftCenter()
        inPSubtomos = protEdit.inReParticles.get()
        outPSubtomos = getattr(protEdit, editStarOutputs.relionParticles.name, None)
        for inPSubtomo, outPSubtomo in zip(inPSubtomos, outPSubtomos):
            isx, isy, isz = self._getShiftsFromPSubtomogram(inPSubtomo)
            osx, osy, osz = self._getShiftsFromPSubtomogram(outPSubtomo)
            self.assertTrue(abs((isx - 4) - osx) < self.editTestsTol)
            self.assertTrue(abs((isy - 2) - osy) < self.editTestsTol)
            self.assertTrue(abs((isz - 3) - osz) < self.editTestsTol)

    def testEditStar_addToAngles(self):
        # Values edited: 5 degrees were added to the rot angle
        addedValue = 5
        protEdit = self._editStar_addToAngles()
        inPSubtomos = protEdit.inReParticles.get()
        outPSubtomos = getattr(protEdit, editStarOutputs.relionParticles.name, None)

        for inPSubtomo, outPSubtomo in zip(inPSubtomos, outPSubtomos):
            irot, itilt, ipsi = self._getAnglesFromPSubtomogram(inPSubtomo)
            orot, otilt, opsi = self._getAnglesFromPSubtomogram(outPSubtomo)
            newRot = irot + addedValue
            # Angle must be expressed in a range of [-180, 180]
            if newRot > 180:
                expectedRot = -180 + (newRot - 180)
            elif newRot < -180:
                expectedRot = 180 + (newRot + 180)
            else:
                expectedRot = newRot
            self.assertTrue(abs(orot - expectedRot) < self.editTestsTol)
            self.assertTrue(abs(itilt - otilt) < self.editTestsTol)
            self.assertTrue(abs(ipsi - opsi) < self.editTestsTol)

    def testEditStar_multiplyCoordinates(self):
        # Values edited: multiply by 2 the X and Z coordinates
        val = 2
        protEdit = self._editStar_multiplyCoordinates()
        inPSubtomos = protEdit.inReParticles.get()
        outPSubtomos = getattr(protEdit, editStarOutputs.relionParticles.name, None)
        for inPSubtomo, outPSubtomo in zip(inPSubtomos, outPSubtomos):
            ix, iy, iz = inPSubtomo.getCoords()
            ox, oy, oz = outPSubtomo.getCoords()
            self.assertTrue(abs(ix * val - ox) < self.editTestsTol)
            self.assertTrue(abs(iy - oy) < self.editTestsTol)
            self.assertTrue(abs(iz * val - oz) < self.editTestsTol)

    def testEditStar_setCoordinatesToValue(self):
        # Values edited: set the Y and Z coordinates to 123
        val = 123
        protEdit = self._editStar_setCoordinatesToValue()
        inPSubtomos = protEdit.inReParticles.get()
        outPSubtomos = getattr(protEdit, editStarOutputs.relionParticles.name, None)
        for inPSubtomo, outPSubtomo in zip(inPSubtomos, outPSubtomos):
            ix, iy, iz = inPSubtomo.getCoords()
            ox, oy, oz = outPSubtomo.getCoords()
            self.assertTrue(abs(ix - ox) < self.editTestsTol)
            self.assertTrue(abs(oy - val) < self.editTestsTol)
            self.assertTrue(abs(oz - val) < self.editTestsTol)

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
    def _editStar_addToAngles(cls) -> Type[ProtRelionEditParticlesStar]:
        print(magentaStr("\n==> Perform angle re-assignment:"))
        protEditStar = cls.newProtocol(ProtRelionEditParticlesStar,
                                       inReParticles=cls.inReParticles,
                                       doRecenter=False,
                                       chosenOperation=cls.editStarOperationDict[OP_ADDITION],
                                       opValue=5,
                                       operateWith=cls.editStarLabelsDict[ANGLES],
                                       label1rot=True)
        protEditStar.setObjLabel('Edit angles')
        cls.launchProtocol(protEditStar)
        return protEditStar

    @classmethod
    def _editStar_multiplyCoordinates(cls) -> Type[ProtRelionEditParticlesStar]:
        print(magentaStr("\n==> Perform coordinates multiply by a scalar:"))
        protEditStar = cls.newProtocol(ProtRelionEditParticlesStar,
                                       inReParticles=cls.inReParticles,
                                       doRecenter=False,
                                       chosenOperation=cls.editStarOperationDict[OP_MULTIPLICATION],
                                       opValue=2,
                                       operateWith=cls.editStarLabelsDict[COORDINATES],
                                       label1x=True,
                                       label3z=True)
        protEditStar.setObjLabel('Multiply coords')
        cls.launchProtocol(protEditStar)
        return protEditStar

    @classmethod
    def _editStar_setCoordinatesToValue(cls) -> Type[ProtRelionEditParticlesStar]:
        print(magentaStr("\n==> Perform coordinates setting to a introduced valuer:"))
        protEditStar = cls.newProtocol(ProtRelionEditParticlesStar,
                                       inReParticles=cls.inReParticles,
                                       doRecenter=False,
                                       chosenOperation=cls.editStarOperationDict[OP_SET_TO],
                                       opValue=123,
                                       operateWith=cls.editStarLabelsDict[COORDINATES],
                                       label2y=True,
                                       label3z=True)
        protEditStar.setObjLabel('Set coords to value')
        cls.launchProtocol(protEditStar)
        return protEditStar

    @classmethod
    def _getShiftsFromPSubtomogram(cls, pSubtomo):
        M = pSubtomo.getTransform().getMatrix()
        sx, sy, sz = translation_from_matrix(M)
        return sx, sy, sz

    @classmethod
    def _getAnglesFromPSubtomogram(cls, pSubtomo):
        angles, _ = getTransformInfoFromCoordOrSubtomo(pSubtomo, pSubtomo.getSamplingRate())
        return angles[:]


class TestRelionTomoGenInitialModel(TestRelion5RefineCycleBase):

    @classmethod
    def _runPreviousProtocols(cls):
        super()._runPreviousProtocols()
        cls.protExtract = cls._runExtractSubtomos()
        cls.inReParticles = getattr(cls.protExtract, cls.protExtract._possibleOutputs.relionParticles.name, None)

    def testInitialModel(self):
        recVol = self._genInitialModel()
        self.checkAverage(recVol,
                          expectedSRate=DataSetRe4STATuto.unbinnedPixSize.value * self.binFactor4,
                          expectedBoxSize=self.croppedBoxSizeBin6,
                          hasHalves=False)

    @classmethod
    def _genInitialModel(cls):
        print(magentaStr("\n==> Generating the a de novo 3D initial model:"))
        protInitialModel = cls.newProtocol(ProtRelionDeNovoInitialModel,
                                           inReParticles=cls.inReParticles,
                                           nVdamMiniBatches=10,
                                           maskDiameter=230,
                                           symmetry=cls.symmetry,
                                           doInC1AndApplySymLater=False,
                                           pooledSubtomos=3,
                                           doGpu=True,
                                           gpusToUse='0',
                                           numberOfMpi=1,
                                           numberOfThreads=3)
        cls.launchProtocol(protInitialModel)
        return getattr(protInitialModel, ProtRelionDeNovoInitialModel._possibleOutputs.average.name, None)


class TestRelionTomo3dClassify(TestRelion5RefineCycleBase):
    nClasses = 2

    @classmethod
    def _runPreviousProtocols(cls):
        super()._runPreviousProtocols()
        cls.importedRef = cls._runImportReference()
        cls.protExtract = cls._runExtractSubtomos()
        cls.inReParticles = getattr(cls.protExtract, cls.protExtract._possibleOutputs.relionParticles.name, None)

    def testCl3d(self):
        protCl3d = self._run3dClassify()
        self._checkCl3dResults(protCl3d)

    def testCl3dWithAlignment(self):
        protCl3d = self._run3dClassify(doAlingment=True)
        self._checkCl3dResults(protCl3d, onlyClassify=False)

    @classmethod
    def _run3dClassify(cls, doAlingment=False):
        paramsDict = {'inReParticles': cls.inReParticles,
                      'referenceVolume': cls.importedRef,
                      'numberOfClasses': cls.nClasses,
                      'initialLowPassFilterA': 30,
                      'symmetry': cls.symmetry,
                      'maskDiameter': 230,
                      'pooledSubtomos': 6,
                      'numberOfMpi': 1,
                      'numberOfThreads': 3}

        if doAlingment:
            paramsDict['doImageAlignment'] = True
            paramsDict['nIterations'] = 1  # Prevent error "No orientation was found as better than any other"
            # because of the test dataset
            paramsDict['doGpu'] = True
            paramsDict['gpusToUse'] = '0'
            label = 'cl3d with alignment'
            print(magentaStr("\n==> Classifying and aligning the psudosubtomograms:"))
        else:
            paramsDict['nIterations'] = 3
            label = 'cl3d'
            print(magentaStr("\n==> Classifying the psudosubtomograms:"))

        protCl3d = cls.newProtocol(ProtRelion3DClassifySubtomograms, **paramsDict)
        protCl3d.setObjLabel(label)
        cls.launchProtocol(protCl3d)
        return protCl3d

    def _checkCl3dResults(self, protCl3d, onlyClassify=True):
        relionParticles = getattr(protCl3d, ProtRelion3DClassifySubtomograms._possibleOutputs.relionParticles.name,
                                  None)
        outClasses = getattr(protCl3d, ProtRelion3DClassifySubtomograms._possibleOutputs.classes.name, None)
        # Check RelionTomoMetadata: only the particles file is generated
        self._checkRe4Metadata(relionParticles,
                               tomogramsFile=self.protExtract._getExtraPath(OUT_TOMOS_STAR),
                               particlesFile=protCl3d._getExtraPath(OUT_PARTICLES_STAR),
                               trajectoriesFile=None,
                               manifoldsFile=None,
                               referenceFscFile=None,
                               relionBinning=self.binFactor4
                               )
        # Check the set of pseudosubtomograms
        expectedSetSize = DataSetRe4STATuto.nCoordsTotal.value
        expectedSRate = DataSetRe4STATuto.unbinnedPixSize.value * self.binFactor4
        expectedBoxSize = self.croppedBoxSizeBin6
        if onlyClassify:
            inCoords = self.protExtract.inputCoords.get()
            self._checkPseudosubtomograms(inCoords, relionParticles,
                                          expectedSetSize=expectedSetSize,
                                          expectedSRate=expectedSRate,
                                          expectedBoxSize=expectedBoxSize,
                                          orientedParticles=True)  # They wer picked with PySeg
        else:
            self.checkRefinedSubtomograms(self.inReParticles, relionParticles,
                                          expectedSetSize=expectedSetSize,
                                          expectedSRate=expectedSRate,
                                          expectedBoxSize=expectedBoxSize,
                                          orientedParticles=True)  # They wer picked with PySeg

        # Check the classes
        self.checkClasses(outClasses,
                          expectedSetSize=self.nClasses,
                          expectedSRate=expectedSRate)


class TestRelionTomoRefine(TestRelion5RefineCycleBase):

    @classmethod
    def _runPreviousProtocols(cls):
        super()._runPreviousProtocols()
        cls.importedRef = cls._runImportReference()
        cls.protExtract = cls._runExtractSubtomos()
        protMakePSubtomos = cls._makePSubtomograms()
        cls.inReParticles = getattr(protMakePSubtomos,
                                    ProtRelionMakePseudoSubtomograms._possibleOutputs.relionParticles.name, None)

    def testAutoRefine(self):
        protAutoRefine = self._runAutoRefine()
        relionParticles = getattr(protAutoRefine, ProtRelionRefineSubtomograms._possibleOutputs.relionParticles.name,
                                  None)
        # Check RelionTomoMetadata: only the particles file is generated
        self._checkRe4Metadata(relionParticles,
                               tomogramsFile=self.protExtract._getExtraPath(OUT_TOMOS_STAR),
                               particlesFile=protAutoRefine._getExtraPath('_data.star'),
                               trajectoriesFile=None,
                               manifoldsFile=None,
                               referenceFscFile=None,
                               relionBinning=self.binFactor4)

        # Check the set of pseudosubtomograms
        expectedSetSize = DataSetRe4STATuto.nCoordsTotal.value
        expectedSRate = DataSetRe4STATuto.unbinnedPixSize.value * self.binFactor4
        expectedBoxSize = self.croppedBoxSizeBin6

        self.checkRefinedSubtomograms(self.inReParticles, relionParticles,
                                      expectedSetSize=expectedSetSize,
                                      expectedSRate=expectedSRate,
                                      expectedBoxSize=expectedBoxSize,
                                      orientedParticles=True)

        # Check the output volume
        avg = getattr(protAutoRefine, ProtRelionRefineSubtomograms._possibleOutputs.average.name, None)
        self.checkAverage(avg,
                          expectedSRate=expectedSRate,
                          expectedBoxSize=expectedBoxSize,
                          hasHalves=True)

    @classmethod
    def _runAutoRefine(cls):
        print(magentaStr("\n==> Refining the particles:"))
        protAutoRefine = cls.newProtocol(ProtRelionRefineSubtomograms,
                                         inReParticles=cls.inReParticles,
                                         referenceVolume=cls.importedRef,
                                         isMapAbsoluteGreyScale=False,
                                         initialLowPassFilterA=50,
                                         symmetry=cls.symmetry,
                                         maskDiameter=230,
                                         useFinerAngularSampling=True,
                                         pooledSubtomos=6,
                                         doGpu=True,
                                         gpusToUse='0',
                                         numberOfMpi=3,
                                         numberOfThreads=3)
        cls.launchProtocol(protAutoRefine)
        return protAutoRefine
