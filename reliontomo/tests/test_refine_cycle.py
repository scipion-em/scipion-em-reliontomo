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

from cistem.protocols import CistemProtTsCtffind
from pyworkflow.tests import BaseTest, setupTestProject, DataSet
from pyworkflow.utils import magentaStr
from reliontomo.constants import OUT_TOMOS_STAR, OUT_PARTICLES_STAR
from reliontomo.protocols import ProtImportCoordinates3DFromStar, ProtRelionPrepareData, \
    ProtRelionMakePseudoSubtomograms, ProtRelionDeNovoInitialModel, ProtRelionRefineSubtomograms, \
    ProtRelionReconstructParticle
from reliontomo.protocols.protocol_make_pseudo_subtomos import outputObjects as makePSubtomosOutputs
from reliontomo.protocols.protocol_prepare_data import outputObjects as prepareOutputs
from reliontomo.protocols.protocol_de_novo_initial_model import outputObjects as iniModelOutputs
from reliontomo.protocols.protocol_prepare_data import ETOMO_FROM_DIR
from reliontomo.tests import RE4_TOMO, DataSetRe4Tomo, OUTPUT_TOMOS, OUTPUT_COORDS
from tomo.protocols import ProtImportTomograms, ProtImportTs

RELION_TOMO_MD = prepareOutputs.outputRelionParticles.name
OUTPUT_VOLUMES = makePSubtomosOutputs.outputVolumes.name
OUTPUT_MODEL = iniModelOutputs.outputAverage.name


class TestRefinceCycle(BaseTest):
    protAutoRefine = None
    protInitialModel = None
    protMakePSubtomos = None
    protPrepare = None
    inCoords = None
    ctfTomoSeries = None
    inTS = None
    inTomoSet = None
    dataset = None
    nParticles = 56
    boxSizeBin4 = 96
    boxSizeBin2 = 128
    samplingRate = 1.35
    tsId = 'TS_43'
    symmetry = 'C6'

    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.dataset = DataSet.getDataSet(RE4_TOMO)
        cls.inTomoSet = cls._importTomograms()
        cls.inCoords = cls._importCoords3dFromStarFile()
        cls.inTS = cls._importTS()
        cls.ctfTomoSeries = cls._estimateCTF()
        cls.protPrepare = cls._prepareData4RelionTomo()
        cls.protMakePSubtomos = cls._makePSubtomograms()
        cls.protInitialModel = cls._genInitialModel()
        cls.protAutoRefine = cls._autoRefine()
        cls.protRecPartFromTS = cls._recParticleFromTS()

    @classmethod
    def _importTomograms(cls):
        print(magentaStr("\n==> Importing data - tomograms:"))
        tomogramsBinning = 4
        protImportTomogram = cls.newProtocol(ProtImportTomograms,
                                             filesPath=cls.dataset.getFile(DataSetRe4Tomo.tomogram.name),
                                             samplingRate=tomogramsBinning * cls.samplingRate)

        cls.launchProtocol(protImportTomogram)
        outputTomos = getattr(protImportTomogram, OUTPUT_TOMOS, None)
        return outputTomos

    @classmethod
    def _importCoords3dFromStarFile(cls):
        print(magentaStr("\n==> Importing data - coordinates from star file:"))
        protImportCoords3dFromStar = cls.newProtocol(ProtImportCoordinates3DFromStar,
                                                     starFile=cls.dataset.getFile(DataSetRe4Tomo.coordinates.name),
                                                     inTomos=cls.inTomoSet,
                                                     samplingRate=1.35,
                                                     boxSize=cls.boxSizeBin4)

        cls.launchProtocol(protImportCoords3dFromStar)
        outCoords = getattr(protImportCoords3dFromStar, OUTPUT_COORDS, None)
        return outCoords

    @classmethod
    def _importTS(cls):
        print(magentaStr("\n==> Importing data - tilt series:"))
        protImportTS = cls.newProtocol(ProtImportTs,
                                       filesPath=cls.dataset.getFile(DataSetRe4Tomo.eTomoDir.name),
                                       filesPattern='*/*.mdoc',
                                       voltage=300,
                                       ampContrast=0.07,
                                       samplingRate=cls.samplingRate,
                                       dosePerFrame=3.05)
        cls.launchProtocol(protImportTS)
        outputTS = getattr(protImportTS, 'outputTiltSeries', None)
        return outputTS

    @classmethod
    def _estimateCTF(cls):
        print(magentaStr("\n==> Estimating the CTF with cistem - ctffind:"))
        protCtfEst = cls.newProtocol(CistemProtTsCtffind,
                                     inputTiltSeries=cls.inTS,
                                     numberOfThreads=6)
        cls.launchProtocol(protCtfEst)
        outputTS = getattr(protCtfEst, 'outputSetOfCTFTomoSeries', None)
        return outputTS

    @classmethod
    def _prepareData4RelionTomo(cls):
        print(magentaStr("\n==> Preparing data for relion 4 tomo:"))
        protPrepare = cls.newProtocol(ProtRelionPrepareData,
                                      inputCtfTs=cls.ctfTomoSeries,
                                      eTomoDataFrom=ETOMO_FROM_DIR,
                                      eTomoFilesPath=cls.dataset.getFile(DataSetRe4Tomo.eTomoDir.name),
                                      flipYZ=True,
                                      flipZ=True,
                                      inputCoords=cls.inCoords)
        cls.launchProtocol(protPrepare)
        return protPrepare

    @classmethod
    def _makePSubtomograms(cls):
        print(magentaStr("\n==> Making the psudosubtomograms:"))
        protMakePsubtomos = cls.newProtocol(ProtRelionMakePseudoSubtomograms,
                                            inOptSet=getattr(cls.protPrepare, RELION_TOMO_MD, None),
                                            boxSize=192,
                                            croppedBoxSize=cls.boxSizeBin4,
                                            binningFactor=4,
                                            outputInFloat16=False,
                                            numberOfThreads=5,
                                            numberOfMpi=3)
        cls.launchProtocol(protMakePsubtomos)
        return protMakePsubtomos

    @classmethod
    def _genInitialModel(cls):
        print(magentaStr("\n==> Generating the a de novo 3D initial model:"))
        protInitialModel = cls.newProtocol(ProtRelionDeNovoInitialModel,
                                           inOptSet=getattr(cls.protMakePSubtomos, RELION_TOMO_MD, None),
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
    def _autoRefine(cls):
        print(magentaStr("\n==> Refining the particles:"))
        protAutoRefine = cls.newProtocol(ProtRelionRefineSubtomograms,
                                         inOptSet=getattr(cls.protMakePSubtomos, RELION_TOMO_MD, None),
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
    def _recParticleFromTS(cls):
        print(magentaStr("\n==> Reconstructing the particle from the TS using a binning factor of 2:"))
        protRecPartFromTS = cls.newProtocol(ProtRelionReconstructParticle,
                                            inOptSet=getattr(cls.protAutoRefine, RELION_TOMO_MD, None),
                                            boxSize=256,
                                            croppedBoxSize=cls.boxSizeBin2,
                                            binningFactor=2,
                                            symmetry=cls.symmetry,
                                            outputInFloat16=False,
                                            numberOfThreads=5,
                                            numberOfMpi=3)
        cls.launchProtocol(protRecPartFromTS)
        return protRecPartFromTS

    def _checkRe4Metadata(self, mdObj, tomogramsFile=None, particlesFile=None, trajectoriesFile=None,
                          manifoldsFile=None, referenceFscFile=None, relionBinning=None):
        self.assertEqual(self.nParticles, mdObj.getNumParticles())
        self.assertEqual(self.samplingRate, mdObj.getTsSamplingRate())
        self.assertEqual(tomogramsFile, mdObj.getTomograms())
        self.assertEqual(particlesFile, mdObj.getParticles())
        self.assertEqual(trajectoriesFile, mdObj.getTrajectories())
        self.assertEqual(manifoldsFile, mdObj.getManifolds())
        self.assertEqual(referenceFscFile, mdObj.getReferenceFsc())
        self.assertEqual(relionBinning, mdObj.getRelionBinning())
        self.assertEqual(self.samplingRate * relionBinning, mdObj.getCurrentSamplingRate())

    def _checkPseudosubtomograms(self, pSubtomosSet, boxSize=None, currentSRate=None):
        self.assertSetSize(pSubtomosSet, self.nParticles)
        self.assertEqual(currentSRate, pSubtomosSet.getSamplingRate())
        for pSubtomo in pSubtomosSet:
            self.assertTrue(exists(pSubtomo.getFileName()))
            self.assertTrue(exists(pSubtomo.getCtfFile()))
            self.assertEqual((boxSize, boxSize, boxSize), pSubtomo.getDimensions())
            self.assertEqual(currentSRate, pSubtomosSet.getSamplingRate())
            self.assertEqual(self.tsId, pSubtomo.getTomoId())

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
        self.assertIsNotNone(self.inTomoSet, 'No tomograms were genetated.')
        self.assertIsNotNone(self.inCoords, 'No coordinates were genetated.')
        self.assertIsNotNone(self.inTS, 'No tilt series were genetated.')
        self.assertIsNotNone(self.ctfTomoSeries, 'No CTF tomo series were genetated.')

    def testPrepareData(self):
        protPrepare = self.protPrepare
        # Check RelionTomoMetadata: both particles and tomograms files are generated
        self._checkRe4Metadata(getattr(protPrepare, RELION_TOMO_MD, None),
                               tomogramsFile=protPrepare._getExtraPath(OUT_TOMOS_STAR),
                               particlesFile=protPrepare._getExtraPath(OUT_PARTICLES_STAR),
                               trajectoriesFile=None,
                               manifoldsFile=None,
                               referenceFscFile=None,
                               relionBinning=4  # Tomograms and coordinates are at bin 4 respecting the TS
                               )

    def testMakePSubtomos(self):
        protMakePSubtomos = self.protMakePSubtomos
        mdObj = getattr(protMakePSubtomos, RELION_TOMO_MD, None)
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
        self._checkPseudosubtomograms(getattr(protMakePSubtomos, OUTPUT_VOLUMES, None),
                                      boxSize=protMakePSubtomos.croppedBoxSize.get(),
                                      currentSRate=mdObj.getCurrentSamplingRate())

    def testInitialModel(self):
        protInitialModel = self.protInitialModel
        recVol = getattr(protInitialModel, OUTPUT_MODEL, None)
        self._checkRecVolume(recVol, optSet=protInitialModel.inOptSet.get(), boxSize=self.boxSizeBin4)

    def testAutoRefine(self):
        protMakePSubtomos = self.protMakePSubtomos
        protAutoRefine = self.protAutoRefine
        mdObj = getattr(protAutoRefine, RELION_TOMO_MD, None)
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
        self._checkPseudosubtomograms(getattr(protAutoRefine, OUTPUT_VOLUMES, None),
                                      boxSize=protMakePSubtomos.croppedBoxSize.get(),
                                      currentSRate=mdObj.getCurrentSamplingRate()
                                      )
        # Check the output volume
        recVol = getattr(protAutoRefine, OUTPUT_MODEL, None)
        pattern = '*it*half%s_class*.mrc'
        half1 = self._getLastFileName(protAutoRefine._getExtraPath(pattern % 1))
        half2 = self._getLastFileName(protAutoRefine._getExtraPath(pattern % 2))
        self._checkRecVolume(recVol,
                             optSet=protAutoRefine.inOptSet.get(),
                             boxSize=self.boxSizeBin4,
                             halves=[half1, half2])

    def testRecParticleFromTS(self):
        protRecPartFromTS = self.protRecPartFromTS
        mdObj = getattr(protRecPartFromTS, RELION_TOMO_MD, None)
        # Check RelionTomoMetadata: only the particles file is generated
        self._checkRe4Metadata(mdObj,
                               tomogramsFile=self.protPrepare._getExtraPath(OUT_TOMOS_STAR),
                               particlesFile=self.protAutoRefine._getExtraPath(OUT_PARTICLES_STAR),
                               trajectoriesFile=None,
                               manifoldsFile=None,
                               referenceFscFile=None,
                               relionBinning=4
                               )
        # Check the output volume
        recVol = getattr(protRecPartFromTS, OUTPUT_MODEL, None)
        self._checkRecVolume(recVol,
                             optSet=protRecPartFromTS.inOptSet.get(),
                             boxSize=self.boxSizeBin2,
                             halves=[protRecPartFromTS._getExtraPath('half1.mrc'),
                                     protRecPartFromTS._getExtraPath('half2.mrc')])



