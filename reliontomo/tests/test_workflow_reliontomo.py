# **************************************************************************
# *
# * Authors:     Scipion Team
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 3 of the License, or
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

import os
from os.path import exists

import pyworkflow.tests as pwtests
from imod.protocols import ProtImodAutomaticCtfEstimation
from pwem.tests.workflows import TestWorkflow
from pwem.protocols import ProtImportVolumes, ProtImportMask
from pyworkflow.utils import magentaStr
from reliontomo.protocols.protocol_ctf_3d_estimation import CTF3D_PER_SUBVOLUME
from tomo.protocols import ProtImportTomograms, ProtImportCoordinates3D, ProtImportTs
from dynamo.protocols import DynamoExtraction

from tomo.protocols.protocol_import_coordinates import IMPORT_FROM_DYNAMO, IMPORT_FROM_CHOICES
from ..protocols import *


CPUS = os.environ.get('SCIPION_TEST_CPUS', 4)
GPUS = os.environ.get('SCIPION_TEST_GPUS', 2)


class TestWorkflowRelionTomo(TestWorkflow):

    boxSize = None
    ds = None
    samplingRate = None
    protImportTomo = None
    protImportCoords3D = None
    protImportTS = None
    protTSCtfEstimationImod = None
    protEstimateCTF3D = None
    protExtractSubtomo = None
    protImportRefVol = None
    protImportMask = None
    protClassify = None
    protRefine = None
    protReconstruct = None

    @classmethod
    def setUpClass(cls):
        pwtests.setupTestProject(cls)
        cls.ds = pwtests.DataSet.getDataSet('reliontomo')
        cls.samplingRate = 4.4
        cls.boxSize = 72
        cls.protImportTomo = cls._importTomograms()
        cls.protImportCoords3D = cls._importCoordinates3D()
        cls.protImportTS = cls._importTiltSeries()
        cls.protTSCtfEstimationImod = cls._tSCtfEstimateImod()
        cls.protEstimateCTF3D = cls._estimateCTF3D()
        cls.protImportRefVol = cls._importVolume()
        cls.protImportMask = cls._importMask()
        cls.protExtractSubtomo = cls._extractSubtomograms()

    @ classmethod
    def _importTomograms(cls):
        print(magentaStr("\n==> Importing data - tomograms:"))
        protImportTomo = cls.newProtocol(
            ProtImportTomograms,
            filesPath=cls.ds.getPath(),
            filesPattern='*Binned1.mrc',
            samplingRate=cls.samplingRate
        )
        protImportTomo.setObjLabel('import 2 tomograms')
        protImportTomo = cls.launchProtocol(protImportTomo)

        return protImportTomo

    @classmethod
    def _importCoordinates3D(cls):
        print(magentaStr("\n==> Importing data - coordinates 3D:"))
        protImportCoords3D = cls.newProtocol(
            ProtImportCoordinates3D,
            importFrom=IMPORT_FROM_CHOICES.index(IMPORT_FROM_DYNAMO),
            filesPath=cls.ds.getPath(),
            filesPattern='*.tbl',
            samplingRate=cls.samplingRate,
            boxSize=cls.boxSize,
            importTomograms=getattr(cls.protImportTomo, 'outputTomograms', None)
        )
        protImportCoords3D.setObjLabel('import 3D coordinates')
        protImportCoords3D = cls.launchProtocol(protImportCoords3D)

        return protImportCoords3D

    @classmethod
    def _importTiltSeries(cls):
        print(magentaStr("\n==> Importing data - tilt series:"))
        protImportTS = cls.newProtocol(
            ProtImportTs,
            filesPath=cls.ds.getPath(),
            filesPattern='{TS}_bin1_ali_bin1.mrcs',
            anglesFrom=2,  # ANGLES_FROM_TLT
            magnification=10000,
            samplingRate=cls.samplingRate,
            voltage=300,
            dosePerFrame=2.1,
            tiltAxisAngle=87.2
        )
        protImportTS.setObjLabel('import tilt series')
        protImportTS = cls.launchProtocol(protImportTS)

        return protImportTS

    @classmethod
    def _tSCtfEstimateImod(cls):
        print(magentaStr("\n==> Calculating the tilt series ctf:"))
        protTSCtfImod = cls.newProtocol(
            ProtImodAutomaticCtfEstimation,
            inputSet=getattr(cls.protImportTS, 'outputTiltSeries', None),
            angleRange=10
        )
        protTSCtfImod.setObjLabel('calculate TS CTF')
        cls.launchProtocol(protTSCtfImod)

        return protTSCtfImod

    @classmethod
    def _estimateCTF3D(cls):
        print(magentaStr("\n==> Estimating the 3D CTF:"))
        protEstimateCTF3D = cls.newProtocol(
            ProtRelionEstimateCTF3D,
            inputCoordinates=getattr(cls.protImportCoords3D, 'outputCoordinates', None),
            inputSetCTFTomoSeries=getattr(cls.protTSCtfEstimationImod, 'outputSetOfCTFTomoSeries', None),
            doseFilesPath=cls.ds.getPath(),
            filesPattern='*ExpDose.txt',
            boxSize=72,
            ctf3dMode=CTF3D_PER_SUBVOLUME,
        )
        protEstimateCTF3D.setObjLabel('Estimate CTF 3D')
        protEstimateCTF3D = cls.launchProtocol(protEstimateCTF3D)

        return protEstimateCTF3D

    def testEstimateCTF3D(self):
        coord3DSet = getattr(self.protEstimateCTF3D, 'outputCoordinates', None)

        # Validate output tomograms
        self.assertSetSize(coord3DSet, size=150)
        self.assertEqual(coord3DSet.getSamplingRate(), self.samplingRate)
        self.assertEqual(coord3DSet.getBoxSize(), self.boxSize)
        # Output coordinates must have an attribute named _3dcftMrcFile, which stores the
        # path of the each ctf3D file
        for coord3d in coord3DSet:
            self.assertTrue(exists(coord3d._3dcftMrcFile.get()))

    @classmethod
    def _importVolume(cls):
        print(magentaStr("\n==> Importing reference volume:"))
        protImportVol = cls.newProtocol(
            ProtImportVolumes,
            emdbId=0,  # IMPORT_FROM_FILES
            filesPath=cls.ds.getFile('refVol'),
            samplingRate=cls.samplingRate,
        )
        protImportVol.setObjLabel('Import reference volume')
        protImportVol = cls.launchProtocol(protImportVol)

        return protImportVol

    @classmethod
    def _importMask(cls):
        print(magentaStr("\n==> Importing the mask:"))
        protImportMask = cls.newProtocol(
            ProtImportMask,
            maskPath=cls.ds.getFile('mask'),
            samplingRate=cls.samplingRate,
        )
        protImportMask.setObjLabel('Import mask')
        protImportMask = cls.launchProtocol(protImportMask)

        return protImportMask

    @classmethod
    def _extractSubtomograms(cls):
        print(magentaStr("\n==> Extracting subtomograms:"))
        protDynamoExtract = cls.newProtocol(
            DynamoExtraction,
            inputCoordinates=getattr(cls.protEstimateCTF3D, 'outputCoordinates', None),
            tomoSource=1,  # OTHER
            boxSize=cls.boxSize,
            doInvert=True,
            inputTomograms=getattr(cls.protImportTomo, 'outputTomograms', None)
        )
        protDynamoExtract.setObjLabel('extract subtomograms')
        protDynamoExtract = cls.launchProtocol(protDynamoExtract)
        return protDynamoExtract

    def testClassifySubtomograms(self):
        print(magentaStr("\n==> Classifying subtomograms:"))
        protClassifSubtomo = self.newProtocol(
            ProtRelionSubtomoClassif3D,
            threads=1,
            mpi=4,
            inputSubtomograms=getattr(self.protExtractSubtomo, 'outputSetOfSubtomogram', None),
            numberOfClasses=1,
            referenceVolume=getattr(self.protImportRefVol, 'outputVolume', None),
            doCTF=False,
            hasReferenceCTFCorrected=True,
            ctfMultiplied=True,
            referenceMask=getattr(self.protImportMask, 'outputMask', None),
            pooledSubtomos=3,
            numberOfIterations=5
        )
        protClassifSubtomo.setObjLabel('classify subtomograms')
        protClassifSubtomo = self.launchProtocol(protClassifSubtomo)
        subTomoClasses = getattr(protClassifSubtomo, 'outputClasses', None)
        avgSubTomoSet = getattr(protClassifSubtomo, 'outputVolumes', None)

        # Validate output classes
        self.assertSetSize(subTomoClasses, size=1)

        # Validate output average subtomograms
        self.assertSetSize(avgSubTomoSet, size=1)
        self.assertEqual(avgSubTomoSet.getSamplingRate(), self.samplingRate)
        self.assertEqual(avgSubTomoSet.getDim(), (self.boxSize, self.boxSize, self.boxSize))

        self.protClassify = protClassifSubtomo

    def testRefineSubtomograms(self):
        print(magentaStr("\n==> Refining subtomograms:"))
        protRefineSubtomo = self.newProtocol(
            ProtRelionSubtomoRefine3D,
            threads=3,
            mpi=5,
            inputSubtomograms=getattr(self.protExtractSubtomo, 'outputSetOfSubtomogram', None),
            numberOfClasses=1,
            referenceVolume=getattr(self.protImportRefVol, 'outputVolume', None),
            doCTF=False,
            hasReferenceCTFCorrected=True,
            ctfMultiplied=True,
            referenceMask=getattr(self.protImportMask, 'outputMask', None),
            pooledSubtomos=10,
            doGpu=True,
            gpusToUse='0',
            numberOfIterations=5,
            extraParams='--sigma_tilt 3.667 --sigma_psi 3.667'
        )
        protRefineSubtomo.setObjLabel('refine subtomogram')
        protRefineSubtomo = self.launchProtocol(protRefineSubtomo)
        subTomoSet = getattr(protRefineSubtomo, 'outputSubtomograms', None)
        avgSubTomo = getattr(protRefineSubtomo, 'outputVolume', None)

        # Validate output subtomograms
        self.assertSetSize(subTomoSet, size=150)
        self.assertEqual(subTomoSet.getSamplingRate(), self.samplingRate)
        self.assertEqual(subTomoSet.getDim(), (self.boxSize, self.boxSize, self.boxSize))

        # Validate output average subtomogram
        self.assertEqual(avgSubTomo.getSamplingRate(), self.samplingRate)
        self.assertEqual(avgSubTomo.getDim(), (self.boxSize, self.boxSize, self.boxSize))

        self.protRefine = protRefineSubtomo

    def testReconstructSubtomograms(self):
        print(magentaStr("\n==> Reconstructing subtomograms:"))
        protReconstructSubtomo = self.newProtocol(
            ProtRelionSubTomoReconstruct,
            threads=3,
            mpi=5,
            inputSubtomos=getattr(self.protRefine, 'outputSetOfSubtomograms', None),
            maxRes=40
        )
        protReconstructSubtomo.setObjLabel('reconstruct subtomograms')
        protReconstructSubtomo = self.launchProtocol(protReconstructSubtomo)
        recTomo = getattr(protReconstructSubtomo, 'outputAvgSubtomogram', None)

        # Validate output average subtomograms
        self.assertEqual(recTomo.getSamplingRate(), self.samplingRate)
        self.assertEqual(recTomo.getDim(), (self.boxSize, self.boxSize, self.boxSize))

    def testCheckTransformationMatrix(self):
        """The output of the refine protocol (convertOutputStep) must be the same as the input of the reconstruct
        protocol (convertInputStep). This way, the bidirectional equivalence can be guaranteed."""

