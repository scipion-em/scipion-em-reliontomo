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
import glob
from os.path import exists, getmtime
from gctf.protocols import ProtTsGctf
from imod.protocols import ProtImodImportTransformationMatrix, ProtImodApplyTransformationMatrix, \
    ProtImodTSNormalization
from pwem.convert.transformations import translation_from_matrix, euler_from_matrix
from pwem.protocols import ProtImportMask
from pwem.protocols.protocol_import.masks import ImportMaskOutput
from pyworkflow.tests import BaseTest, setupTestProject, DataSet
from pyworkflow.utils import magentaStr
from reliontomo.constants import OUT_TOMOS_STAR, OUT_PARTICLES_STAR, IN_PARTICLES_STAR, POSTPROCESS_DIR, \
    POST_PROCESS_MRC
from reliontomo.protocols import ProtImportCoordinates3DFromStar, ProtRelionPrepareData, \
    ProtRelionMakePseudoSubtomograms, ProtRelionDeNovoInitialModel, ProtRelionRefineSubtomograms, \
    ProtRelionReconstructParticle, ProtExtractCoordsFromPSubtomos, ProtRelionTomoReconstruct, \
    ProtRelionEditParticlesStar
from reliontomo.protocols.protocol_3d_classify_subtomograms import outputObjects as cl3dOutputs, \
    ProtRelion3DClassifySubtomograms
from reliontomo.protocols.protocol_base_import_from_star import outputObjects as importStarOutputs
from reliontomo.protocols.protocol_edit_particles_star import OPERATION_LABELS, LABELS_TO_OPERATE_WITH, ANGLES, \
    OP_ADDITION, OP_MULTIPLICATION, COORDINATES, OP_SET_TO
from reliontomo.protocols.protocol_extract_coordinates_from_psubtomos import outputObjects as extractCoordsOutputs
from reliontomo.protocols.protocol_prepare_data import outputObjects as prepareOutputs
from reliontomo.protocols.protocol_edit_particles_star import outputObjects as editStarOutputs
from reliontomo.protocols.protocol_de_novo_initial_model import outputObjects as iniModelOutputs
from reliontomo.protocols.protocol_rec_tomogram import outputObjects as recTomoRelionOutputs
from reliontomo.protocols.protocol_reconstruc_particle_from_ts import outputObjects as recParticleFromTsOutputs
from reliontomo.protocols.protocol_rec_tomogram import SINGLE_TOMO, ALL_TOMOS
from reliontomo.tests import RE4_TOMO, DataSetRe4Tomo
from reliontomo.utils import genEnumParamDict
from tomo.protocols import ProtImportTs
from tomo3d.protocols import ProtJjsoftReconstructTomogram
from tomo3d.protocols.protocol_base_reconstruct import outputTomoRecObjects
from tomo3d.protocols.protocol_reconstruct_tomogram import SIRT

RELION_TOMO_PARTICLES = prepareOutputs.relionParticles.name
OUTPUT_MODEL = iniModelOutputs.average.name
OUTPUT_COORDS = importStarOutputs.coordinates.name
OUTPUT_CLASSES = cl3dOutputs.classes.name


class TestRefineCycle(BaseTest):
    fscMaskBin2 = None
    normTS = None
    alignedTS = None
    tsWithAlignment = None
    protAutoRefine = None
    protInitialModel = None
    protMakePSubtomos = None
    protExtractCoords = None
    protPrepare = None
    inCoords = None
    ctfTomoSeries = None
    inTS = None
    inTomoSet = None
    dataset = None
    nParticles = 108
    boxSizeBin4 = 96
    boxSizeBin2 = 128
    samplingRateOrig = 2.7
    tsIds = ['TS_45', 'TS_54']
    symmetry = 'C6'
    nClasses = 2
    editStarOperationDict = genEnumParamDict(OPERATION_LABELS)
    editStarLabelsDict = genEnumParamDict(LABELS_TO_OPERATE_WITH)
    editTestsTol = 0.01

    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.dataset = DataSet.getDataSet(RE4_TOMO)
        cls.fscMaskBin2 = cls._importMask()
        cls.inTS = cls._importTS()
        cls.ctfTomoSeries = cls._estimateCTF()
        cls.tsWithAlignment = cls._importTransformationMatrix()
        cls.alignedTS = cls._applyTransformationMatrix()
        cls.normTS = cls._noramlizeTS()
        cls.inTomoSet = cls._reconstructTomograms()
        cls.inCoords = cls._importCoords3dFromStarFile()
        cls.protPrepare = cls._prepareData4RelionTomo()
        cls.recTomoFromPrepareSingle = cls._runRecFromPrepare_SingleTomo()
        cls.recTomoFromPrepareAll = cls._runRecFromPrepare_AllTomos()
        cls.protMakePSubtomos = cls._makePSubtomograms()
        cls.protEditStarCenter = cls._editStar_shiftCenter()
        cls.protEditStarAngles = cls._editStar_addToAngles()
        cls.protEditStarCoordsMult = cls._editStar_multiplyCoordinates()
        cls.protEditStarSetCoords = cls._editStar_setCoordinatesToValue()
        # cls.protExtractCoords = cls._extractCoordsFromPSubtomos()
        cls.protInitialModel = cls._genInitialModel()
        cls.protCl3d = cls._3dClassify()
        cls.protCl3dWithAlign = cls._3dClassify(doAlingment=True)
        cls.protAutoRefine = cls._autoRefine()
        cls.protRecPartFromTS = cls._recParticleFromTS()
        cls.protRecPartFromTsWithSolvent = cls._recParticleFromTsWithSolventMask()

    @classmethod
    def _importTS(cls):
        print(magentaStr("\n==> Importing the tilt series:"))
        protImportTS = cls.newProtocol(ProtImportTs,
                                       filesPath=cls.dataset.getFile(DataSetRe4Tomo.eTomoDir.name),
                                       filesPattern=DataSetRe4Tomo.mdocs.value,
                                       voltage=300,
                                       ampContrast=0.07,
                                       samplingRate=cls.samplingRateOrig,
                                       dosePerFrame=3.05)
        cls.launchProtocol(protImportTS)
        outputTS = getattr(protImportTS, 'outputTiltSeries', None)
        return outputTS

    @classmethod
    def _importMask(cls):
        print(magentaStr("\n==> Importing the FSC mask:"))
        protImportMask = cls.newProtocol(ProtImportMask,
                                         maskPath=cls.dataset.getFile(DataSetRe4Tomo.maskFscBin2.name),
                                         samplingRate=cls.samplingRateOrig)
        protImportMask.setObjLabel('Import FSC mask')
        cls.launchProtocol(protImportMask)
        outputMask = getattr(protImportMask, ImportMaskOutput.outputMask.name, None)
        return outputMask

    @classmethod
    def _importTransformationMatrix(cls):
        print(magentaStr("\n==> Importing the transformation matrices:"))
        protImportTransMatrix = cls.newProtocol(ProtImodImportTransformationMatrix,
                                                inputSetOfTiltSeries=cls.inTS,
                                                filesPath=cls.dataset.getFile(DataSetRe4Tomo.eTomoDir.name),
                                                filesPattern=DataSetRe4Tomo.alignments.value)
        cls.launchProtocol(protImportTransMatrix)
        return getattr(protImportTransMatrix, 'TiltSeries', None)

    @classmethod
    def _applyTransformationMatrix(cls):
        print(magentaStr("\n==> Applying the transformation matrices:"))
        protApplyTransMatrix = cls.newProtocol(ProtImodApplyTransformationMatrix,
                                               inputSetOfTiltSeries=cls.tsWithAlignment,
                                               binning=1)
        cls.launchProtocol(protApplyTransMatrix)
        return getattr(protApplyTransMatrix, 'InterpolatedTiltSeries', None)

    @classmethod
    def _noramlizeTS(cls):
        print(magentaStr("\n==> Binning the tilt series to bin 4:"))
        protNormTS = cls.newProtocol(ProtImodTSNormalization,
                                     inputSetOfTiltSeries=cls.alignedTS,
                                     binning=4)
        cls.launchProtocol(protNormTS)
        return getattr(protNormTS, 'TiltSeries', None)

    @classmethod
    def _reconstructTomograms(cls):
        print(magentaStr("\n==> Reconstructing the tomograms:"))
        protRecTomograms = cls.newProtocol(ProtJjsoftReconstructTomogram,
                                           inputSetOfTiltSeries=cls.normTS,
                                           method=SIRT,
                                           height=300)  # Thickness at bin 4
        protRecTomograms = cls.launchProtocol(protRecTomograms, wait=True)
        return getattr(protRecTomograms, outputTomoRecObjects.tomograms.name, None)

    @classmethod
    def _importCoords3dFromStarFile(cls):
        print(magentaStr("\n==> Importing data - coordinates from star file:"))
        protImportCoords3dFromStar = cls.newProtocol(ProtImportCoordinates3DFromStar,
                                                     starFile=cls.dataset.getFile(DataSetRe4Tomo.coordinates.name),
                                                     inTomos=cls.inTomoSet,
                                                     samplingRate=cls.samplingRateOrig,
                                                     boxSize=cls.boxSizeBin4)

        cls.launchProtocol(protImportCoords3dFromStar)
        outCoords = getattr(protImportCoords3dFromStar, OUTPUT_COORDS, None)
        return outCoords

    @classmethod
    def _estimateCTF(cls):
        print(magentaStr("\n==> Estimating the CTF with gctf:"))
        protCtfEst = cls.newProtocol(ProtTsGctf,
                                     inputTiltSeries=cls.inTS,
                                     numberOfThreads=6)
        cls.launchProtocol(protCtfEst)
        return getattr(protCtfEst, 'outputSetOfCTFTomoSeries', None)

    @classmethod
    def _prepareData4RelionTomo(cls):
        print(magentaStr("\n==> Preparing data for relion 4 tomo:"))
        protPrepare = cls.newProtocol(ProtRelionPrepareData,
                                      inputCtfTs=cls.ctfTomoSeries,
                                      inputCoords=cls.inCoords,
                                      flipYZ=True,
                                      flipZ=True)
        return cls.launchProtocol(protPrepare)

    @classmethod
    def _recFromPrepare(cls, recMode, tomoId=None):
        protRelionRec = cls.newProtocol(ProtRelionTomoReconstruct,
                                        protPrepare=cls.protPrepare,
                                        recTomoMode=recMode,
                                        tomoId=tomoId,
                                        binFactor=8)
        cls.launchProtocol(protRelionRec)
        return protRelionRec

    @classmethod
    def _runRecFromPrepare_SingleTomo(cls):
        print(magentaStr("\n==> Reconstructing one of the tomograms with Relion:"))
        protRelionRec = cls._recFromPrepare(SINGLE_TOMO, cls.tsIds[0])  # TS_54
        return protRelionRec

    @classmethod
    def _runRecFromPrepare_AllTomos(cls):
        print(magentaStr("\n==> Reconstructing all the tomograms with Relion:"))
        protRelionRec = cls._recFromPrepare(ALL_TOMOS)
        return protRelionRec

    @classmethod
    def _makePSubtomograms(cls):
        print(magentaStr("\n==> Making the psudosubtomograms:"))
        protMakePsubtomos = cls.newProtocol(ProtRelionMakePseudoSubtomograms,
                                            inReParticles=getattr(cls.protPrepare, RELION_TOMO_PARTICLES, None),
                                            boxSize=192,
                                            croppedBoxSize=cls.boxSizeBin4,
                                            binningFactor=4,
                                            outputInFloat16=False,
                                            numberOfThreads=5,
                                            numberOfMpi=3)
        cls.launchProtocol(protMakePsubtomos)
        return protMakePsubtomos

    @classmethod
    def _editStar_shiftCenter(cls) -> ProtRelionEditParticlesStar:
        print(magentaStr("\n==> Perform centering of particles:"))
        protEditStar = cls.newProtocol(ProtRelionEditParticlesStar,
                                       inReParticles=getattr(cls.protMakePSubtomos, RELION_TOMO_PARTICLES, None),
                                       doRecenter=True,
                                       shiftX=4,
                                       shiftY=2,
                                       shiftZ=3)
        protEditStar.setObjLabel('Particles re-center')
        cls.launchProtocol(protEditStar)
        return protEditStar

    @classmethod
    def _editStar_addToAngles(cls) -> ProtRelionEditParticlesStar:
        print(magentaStr("\n==> Perform angle re-assignment:"))
        protEditStar = cls.newProtocol(ProtRelionEditParticlesStar,
                                       inReParticles=getattr(cls.protMakePSubtomos, RELION_TOMO_PARTICLES, None),
                                       doRecenter=False,
                                       chosenOperation=cls.editStarOperationDict[OP_ADDITION],
                                       opValue=5,
                                       operateWith=cls.editStarLabelsDict[ANGLES],
                                       label1rot=True)
        protEditStar.setObjLabel('Edit angles')
        cls.launchProtocol(protEditStar)
        return protEditStar

    @classmethod
    def _editStar_multiplyCoordinates(cls) -> ProtRelionEditParticlesStar:
        print(magentaStr("\n==> Perform coordinates multiply by a scalar:"))
        protEditStar = cls.newProtocol(ProtRelionEditParticlesStar,
                                       inReParticles=getattr(cls.protMakePSubtomos, RELION_TOMO_PARTICLES, None),
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
    def _editStar_setCoordinatesToValue(cls) -> ProtRelionEditParticlesStar:
        print(magentaStr("\n==> Perform coordinates setting to a introduced valuer:"))
        protEditStar = cls.newProtocol(ProtRelionEditParticlesStar,
                                       inReParticles=getattr(cls.protMakePSubtomos, RELION_TOMO_PARTICLES, None),
                                       doRecenter=False,
                                       chosenOperation=cls.editStarOperationDict[OP_SET_TO],
                                       opValue=123,
                                       operateWith=cls.editStarLabelsDict[COORDINATES],
                                       label2y=True,
                                       label3z=True)
        protEditStar.setObjLabel('Set coords to value')
        cls.launchProtocol(protEditStar)
        return protEditStar

    # @classmethod
    # def _extractCoordsFromPSubtomos(cls):
    #     print(magentaStr("\n==> Generating the a de novo 3D initial model:"))
    #     protExtractCoords = cls.newProtocol(ProtExtractCoordsFromPSubtomos,
    #                                         inReParticles=getattr(cls.protMakePSubtomos, RELION_TOMO_PARTICLES, None))
    #     cls.launchProtocol(protExtractCoords)
    #     return protExtractCoords

    @classmethod
    def _genInitialModel(cls):
        print(magentaStr("\n==> Generating the a de novo 3D initial model:"))
        protInitialModel = cls.newProtocol(ProtRelionDeNovoInitialModel,
                                           inReParticles=getattr(cls.protMakePSubtomos, RELION_TOMO_PARTICLES, None),
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
        return protInitialModel

    @classmethod
    def _3dClassify(cls, doAlingment=False):
        paramsDict = {'inReParticles': getattr(cls.protMakePSubtomos, RELION_TOMO_PARTICLES, None),
                      'referenceVolume': getattr(cls.protInitialModel, OUTPUT_MODEL, None),
                      'numberOfClasses': cls.nClasses,
                      'initialLowPassFilterA': 30,
                      'symmetry': cls.symmetry,
                      'maskDiameter': 230,
                      'nIterations': 3,
                      'pooledSubtomos': 6,
                      'numberOfMpi': 1,
                      'numberOfThreads': 3}

        if doAlingment:
            paramsDict['doImageAlignment'] = True
            paramsDict['doGpu'] = True
            paramsDict['gpusToUse'] = '0'
            label = 'cl3d with alignment'
            print(magentaStr("\n==> Classifying the psudosubtomograms with alignment:"))
        else:
            label = 'cl3d'
            print(magentaStr("\n==> Classifying the psudosubtomograms:"))

        protCl3d = cls.newProtocol(ProtRelion3DClassifySubtomograms, **paramsDict)
        protCl3d.setObjLabel(label)
        cls.launchProtocol(protCl3d)
        return protCl3d

    @classmethod
    def _autoRefine(cls):
        print(magentaStr("\n==> Refining the particles:"))
        protAutoRefine = cls.newProtocol(ProtRelionRefineSubtomograms,
                                         inReParticles=getattr(cls.protMakePSubtomos, RELION_TOMO_PARTICLES, None),
                                         referenceVolume=getattr(cls.protInitialModel, OUTPUT_MODEL, None),
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

    @classmethod
    def _genRecPartFromTsDict(cls):
        return {'inReParticles': getattr(cls.protAutoRefine, RELION_TOMO_PARTICLES, None),
                'boxSize': 256,
                'croppedBoxSize': cls.boxSizeBin2,
                'binningFactor': 2,
                'symmetry': cls.symmetry,
                'outputInFloat16': False,
                'numberOfThreads': 5,
                'numberOfMpi': 3}

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
        paramsDict = cls._genRecPartFromTsDict()
        paramsDict['solventMask'] = cls.fscMaskBin2
        protRecPartFromTS = cls.newProtocol(ProtRelionReconstructParticle, **paramsDict)
        protRecPartFromTS.setObjLabel('rec part from ts with solvent mask')
        cls.launchProtocol(protRecPartFromTS)
        return protRecPartFromTS

    def _checkRe4Metadata(self, pSubtomoSet, tomogramsFile=None, particlesFile=None, trajectoriesFile=None,
                          manifoldsFile=None, referenceFscFile=None, relionBinning=None):
        self.assertEqual(self.nParticles, pSubtomoSet.getNReParticles())
        self.assertEqual(self.samplingRateOrig, pSubtomoSet.getTsSamplingRate())
        self.assertEqual(tomogramsFile, pSubtomoSet.getTomograms())
        self.assertEqual(particlesFile, pSubtomoSet.getParticles())
        self.assertEqual(trajectoriesFile, pSubtomoSet.getTrajectories())
        self.assertEqual(manifoldsFile, pSubtomoSet.getManifolds())
        self.assertEqual(referenceFscFile, pSubtomoSet.getReferenceFsc())
        self.assertEqual(relionBinning, pSubtomoSet.getRelionBinning())
        self.assertEqual(self.samplingRateOrig * relionBinning, pSubtomoSet.getCurrentSamplingRate())

    def _checkPseudosubtomograms(self, pSubtomosSet, boxSize=None, currentSRate=None):
        self.assertSetSize(pSubtomosSet, self.nParticles)
        self.assertEqual(currentSRate, pSubtomosSet.getSamplingRate())
        for pSubtomo in pSubtomosSet:
            self.assertTrue(exists(pSubtomo.getFileName()))
            self.assertTrue(exists(pSubtomo.getCtfFile()))
            self.assertEqual((boxSize, boxSize, boxSize), pSubtomo.getDimensions())
            self.assertEqual(currentSRate, pSubtomosSet.getSamplingRate())
            self.assertTrue(pSubtomo.getTsId() in self.tsIds)

    def _check3dClasses(self, classes, currentSRate=None):
        self.assertSetSize(classes, self.nClasses)
        self.assertEqual(currentSRate, classes.getSamplingRate())

    def _checkRecVolume(self, recVol, optSet=None, boxSize=None, halves=None):
        self.assertEqual(recVol.getSamplingRate(), optSet.getCurrentSamplingRate())
        self.assertEqual(recVol.getDim(), (boxSize, boxSize, boxSize))
        if halves:
            half1, half2 = recVol.getHalfMaps().split(',')
            self.assertEqual(halves, [half1, half2])

    @staticmethod
    def _getLastFileName(pattern):
        files = glob.glob(pattern)
        files.sort(key=getmtime)
        return files[-1]

    def testPrevRequiredDataGeneration(self):
        self.assertIsNotNone(self.inTS, 'No tilt series were generated after importing.')
        self.assertIsNotNone(self.tsWithAlignment, 'No tilt series were generated after importing the transformation '
                                                   'matrix.')
        self.assertIsNotNone(self.alignedTS, 'No tilt series were generated after applying the transformation.')
        self.assertIsNotNone(self.normTS, 'No tilt series were generated after normalization.')
        self.assertIsNotNone(self.inTomoSet, 'No tomograms were reconstructed.')
        self.assertIsNotNone(self.inCoords, 'No coordinates were generated.')
        self.assertIsNotNone(self.ctfTomoSeries, 'No CTF tomo series were generated.')

    def testPrepareData(self):
        protPrepare = self.protPrepare
        # Check RelionTomoMetadata: both particles and tomograms files are generated
        self._checkRe4Metadata(getattr(protPrepare, RELION_TOMO_PARTICLES, None),
                               tomogramsFile=protPrepare._getExtraPath(OUT_TOMOS_STAR),
                               particlesFile=protPrepare._getExtraPath(OUT_PARTICLES_STAR),
                               trajectoriesFile=None,
                               manifoldsFile=None,
                               referenceFscFile=None,
                               relionBinning=1
                               )

    def testRecSingleTomoFromPrep(self):
        expectedSize = 1
        self._checkRelionRecTomos(self.recTomoFromPrepareSingle, expectedSize)

    def testRecAllTomosFromPrep(self):
        expectedSize = 2
        self._checkRelionRecTomos(self.recTomoFromPrepareAll, expectedSize)

    def _checkRelionRecTomos(self, protRec, expectedSize):
        outTomos = getattr(protRec, recTomoRelionOutputs.tomograms.name, None)
        self.assertSetSize(outTomos, expectedSize)
        self.assertEqual(outTomos.getFirstItem().getDimensions(), (464, 480, 150))
        self.assertEqual(outTomos.getSamplingRate(), self.samplingRateOrig * protRec.binFactor.get())

    def testMakePSubtomos(self):
        protMakePSubtomos = self.protMakePSubtomos
        mdObj = getattr(protMakePSubtomos, RELION_TOMO_PARTICLES, None)
        # Check RelionTomoMetadata: only the particles file is generated
        self._checkRe4Metadata(mdObj,
                               tomogramsFile=self.protPrepare._getExtraPath(OUT_TOMOS_STAR),
                               particlesFile=protMakePSubtomos._getExtraPath(OUT_PARTICLES_STAR),
                               trajectoriesFile=None,
                               manifoldsFile=None,
                               referenceFscFile=None,
                               relionBinning=4
                               )
        # Check the set of pseudosubtomograms
        self._checkPseudosubtomograms(getattr(protMakePSubtomos, RELION_TOMO_PARTICLES, None),
                                      boxSize=protMakePSubtomos.croppedBoxSize.get(),
                                      currentSRate=mdObj.getCurrentSamplingRate())

    def testEditStar_shiftCenter(self):
        # Values edited: shiftX = 4, shiftY = 2, shiftZ = 3
        protEdit = self.protEditStarCenter
        inPSubtomos = protEdit.inReParticles.get()
        outPSubtomos = getattr(protEdit, editStarOutputs.relionParticles.name, None)
        for inPSubtomo, outPSubtomo in zip (inPSubtomos, outPSubtomos):
            isx, isy, isz = self._getShiftsFromPSubtomogram(inPSubtomo)
            osx, osy, osz = self._getShiftsFromPSubtomogram(outPSubtomo)
            self.assertTrue(abs((isx - 4) - osx) < self.editTestsTol)
            self.assertTrue(abs((isy - 2) - osy) < self.editTestsTol)
            self.assertTrue(abs((isz - 3) - osz) < self.editTestsTol)

    def testEditStar_addToAngles(self):
        # Values edited: 5 degrees were added to the rot angle
        protEdit = self.protEditStarAngles
        inPSubtomos = protEdit.inReParticles.get()
        outPSubtomos = getattr(protEdit, editStarOutputs.relionParticles.name, None)
        for inPSubtomo, outPSubtomo in zip(inPSubtomos, outPSubtomos):
            irot, itilt, ipsi = self._getAnglesFromPSubtomogram(inPSubtomo)
            orot, otilt, opsi = self._getAnglesFromPSubtomogram(outPSubtomo)
            self.assertTrue(abs(irot + 5 - orot) < self.editTestsTol)
            self.assertTrue(abs(itilt - otilt) < self.editTestsTol)
            self.assertTrue(abs(ipsi - opsi) < self.editTestsTol)

    def testEditStar_multiplyCoordinates(self):
        # Values edited: multiply by 2 the X and Z coordinates
        val = 2
        protEdit = self.protEditStarCoordsMult
        inPSubtomos = protEdit.inReParticles.get()
        outPSubtomos = getattr(protEdit, editStarOutputs.relionParticles.name, None)
        for inPSubtomo, outPSubtomo in zip (inPSubtomos, outPSubtomos):
            ix, iy, iz = inPSubtomo.getCoords()
            ox, oy, oz = outPSubtomo.getCoords()
            self.assertTrue(abs(ix * val - ox) < self.editTestsTol)
            self.assertTrue(abs(iy - oy) < self.editTestsTol)
            self.assertTrue(abs(iz * val - oz) < self.editTestsTol)

    def testEditStar_setCoordinatesToValue(self):
        # Values edited: set the Y and Z coordinates to 123
        val = 123
        protEdit = self.protEditStarSetCoords
        inPSubtomos = protEdit.inReParticles.get()
        outPSubtomos = getattr(protEdit, editStarOutputs.relionParticles.name, None)
        for inPSubtomo, outPSubtomo in zip (inPSubtomos, outPSubtomos):
            ix, iy, iz = inPSubtomo.getCoords()
            ox, oy, oz = outPSubtomo.getCoords()
            self.assertTrue(abs(ix - ox) < self.editTestsTol)
            self.assertTrue(abs(oy - val) < self.editTestsTol)
            self.assertTrue(abs(oz - val) < self.editTestsTol)

    @classmethod
    def _getShiftsFromPSubtomogram(cls, pSubtomo):
        M = pSubtomo.getTransform().getMatrix()
        sx, sy, sz = translation_from_matrix(M)
        return sx, sy, sz

    @classmethod
    def _getAnglesFromPSubtomogram(cls, pSubtomo):
        M = pSubtomo.getTransform().getMatrix()
        rot, tilt, psi = euler_from_matrix(M)
        return rot, tilt, psi

    # def testExtractCoordsFromPSubtomos(self):
    #     protExtractCoords = self.protExtractCoords
    #     outCoords = getattr(protExtractCoords, extractCoordsOutputs.coordinates.name, None)
    #     self.assertEqual(outCoords.getSamplingRate(), 5.4)
    #     self.assertEqual(outCoords.getBoxSize(), self.boxSizeBin4)

    def testInitialModel(self):
        protInitialModel = self.protInitialModel
        recVol = getattr(protInitialModel, OUTPUT_MODEL, None)
        self._checkRecVolume(recVol,
                             optSet=protInitialModel.inReParticles.get(),
                             boxSize=self.boxSizeBin4)

    def testCl3d(self):
        self._runTestCl3d(self.protCl3d)

    def testCl3dWithAlignment(self):
        self._runTestCl3d(self.protCl3dWithAlign)

    def _runTestCl3d(self, protCl3d):
        protMakePSubtomos = self.protMakePSubtomos
        mdObj = getattr(protCl3d, RELION_TOMO_PARTICLES, None)
        # Check RelionTomoMetadata: only the particles file is generated
        self._checkRe4Metadata(mdObj,
                               tomogramsFile=self.protPrepare._getExtraPath(OUT_TOMOS_STAR),
                               particlesFile=protCl3d._getExtraPath(OUT_PARTICLES_STAR),
                               trajectoriesFile=None,
                               manifoldsFile=None,
                               referenceFscFile=None,
                               relionBinning=4
                               )
        # Check the set of pseudosubtomograms
        self._checkPseudosubtomograms(getattr(protCl3d, RELION_TOMO_PARTICLES, None),
                                      boxSize=protMakePSubtomos.croppedBoxSize.get(),
                                      currentSRate=mdObj.getCurrentSamplingRate()
                                      )
        # Check the classes
        self._check3dClasses(getattr(protCl3d, OUTPUT_CLASSES, None),
                             currentSRate=mdObj.getCurrentSamplingRate())

    def testAutoRefine(self):
        protMakePSubtomos = self.protMakePSubtomos
        protAutoRefine = self.protAutoRefine
        mdObj = getattr(protAutoRefine, RELION_TOMO_PARTICLES, None)
        # Check RelionTomoMetadata: only the particles file is generated
        self._checkRe4Metadata(mdObj,
                               tomogramsFile=self.protPrepare._getExtraPath(OUT_TOMOS_STAR),
                               particlesFile=protAutoRefine._getExtraPath(OUT_PARTICLES_STAR),
                               trajectoriesFile=None,
                               manifoldsFile=None,
                               referenceFscFile=None,
                               relionBinning=4
                               )
        # Check the set of pseudosubtomograms
        self._checkPseudosubtomograms(getattr(protAutoRefine, RELION_TOMO_PARTICLES, None),
                                      boxSize=protMakePSubtomos.croppedBoxSize.get(),
                                      currentSRate=mdObj.getCurrentSamplingRate()
                                      )
        # Check the output volume
        pattern = '*it*half%s_class*.mrc'
        half1 = self._getLastFileName(protAutoRefine._getExtraPath(pattern % 1))
        half2 = self._getLastFileName(protAutoRefine._getExtraPath(pattern % 2))
        self._checkRecVolume(getattr(protAutoRefine, OUTPUT_MODEL, None),
                             optSet=getattr(protAutoRefine, RELION_TOMO_PARTICLES, None),
                             boxSize=self.boxSizeBin4,
                             halves=[half1, half2])

    def _checkRecVolFromTs(self, protRecPartFromTS, recVol):
        self._checkRecVolume(recVol,
                             optSet=getattr(protRecPartFromTS, RELION_TOMO_PARTICLES, None),
                             boxSize=self.boxSizeBin2,
                             halves=[protRecPartFromTS._getExtraPath('half1.mrc'),
                                     protRecPartFromTS._getExtraPath('half2.mrc')])

    def _checkPSubtomosFromTS(self, protRecPartFromTS, fscFile=None):
        self._checkRe4Metadata(getattr(protRecPartFromTS, RELION_TOMO_PARTICLES, None),
                               tomogramsFile=self.protPrepare._getExtraPath(OUT_TOMOS_STAR),
                               particlesFile=protRecPartFromTS._getExtraPath(IN_PARTICLES_STAR),
                               trajectoriesFile=None,
                               manifoldsFile=None,
                               referenceFscFile=fscFile,
                               relionBinning=2
                               )

    def testRecParticleFromTS(self):
        protRecPartFromTS = self.protRecPartFromTS
        # Check the generated volume
        self._checkRecVolFromTs(protRecPartFromTS, getattr(protRecPartFromTS, OUTPUT_MODEL, None))
        # Check the generated psubtomos
        self._checkPSubtomosFromTS(protRecPartFromTS)

    def testRecParticleFromTsWithSolventMask(self):
        protRecPartFromTsWithSolvent = self.protRecPartFromTsWithSolvent
        # Check the generated volume
        recVol = getattr(protRecPartFromTsWithSolvent, OUTPUT_MODEL, None)
        self._checkRecVolFromTs(protRecPartFromTsWithSolvent, recVol)
        # Check the generated psubtomos: when a solvent mask is introduced, a fsc mask is also generated
        fscFile = protRecPartFromTsWithSolvent._getExtraPath(POSTPROCESS_DIR, POST_PROCESS_MRC.replace('.mrc', '.star'))
        self._checkPSubtomosFromTS(protRecPartFromTsWithSolvent, fscFile=fscFile)
        # Check the generated FSC mask
        fscSolvent = getattr(protRecPartFromTsWithSolvent, recParticleFromTsOutputs.postProcessVolume.name, None)
        self._checkRecVolFromTs(protRecPartFromTsWithSolvent, fscSolvent)

