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
from ..protocols import *


CPUS = os.environ.get('SCIPION_TEST_CPUS', 4)
GPUS = os.environ.get('SCIPION_TEST_GPUS', 2)


class TestWorkflowRelionTomo(TestWorkflow):
    @classmethod
    def setUpClass(cls):
        pwtests.setupTestProject(cls)
        cls.ds = pwtests.DataSet.getDataSet('reliontomo')
        cls.samplingRate = 4.4
        cls.boxSize = 72

    def _importTomograms(self):
        print(magentaStr("\n==> Importing data - tomograms:"))
        protImportTomo = self.newProtocol(
            ProtImportTomograms,
            filesPath=self.ds.getPath(),
            filesPattern='*Binned1.mrc',
            samplingRate=self.samplingRate
        )
        protImportTomo.setObjLabel('import 2 tomograms')
        protImportTomo = self.launchProtocol(protImportTomo)
        tomoSet = getattr(protImportTomo, 'outputTomograms', None)

        # Validate output tomograms
        self.assertSetSize(tomoSet, size=2)
        self.assertEqual(tomoSet.getSamplingRate(), self.samplingRate)
        self.assertEqual(tomoSet.getDim(), (1919, 1854, 500))

        return protImportTomo

    def _importCoordinates3D(self, protImportTomo):
        print(magentaStr("\n==> Importing data - coordinates 3D:"))
        protImportCoords3D = self.newProtocol(
            ProtImportCoordinates3D,
            importFrom=3,  # IMPORT_FROM_DYNAMO
            filesPath=self.ds.getPath(),
            filesPattern='*.tbl',
            samplingRate=self.samplingRate,
            boxSize=self.boxSize,
            importTomograms=getattr(protImportTomo, 'outputTomograms', None)
        )
        protImportCoords3D.setObjLabel('import 3D coordinates')
        protImportCoords3D = self.launchProtocol(protImportCoords3D)
        coord3DSet = getattr(protImportCoords3D, 'outputCoordinates', None)

        # Validate output tomograms
        self.assertSetSize(coord3DSet, size=150)
        self.assertEqual(coord3DSet.getSamplingRate(), self.samplingRate)
        self.assertEqual(coord3DSet.getBoxSize(), self.boxSize)

        return protImportCoords3D

    def _importTiltSeries(self):
        print(magentaStr("\n==> Importing data - tilt series:"))
        protImportTS = self.newProtocol(
            ProtImportTs,
            filesPath=self.ds.getPath(),
            filesPattern='{TS}_bin1_ali_bin1.mrcs',
            anglesFrom=2,  # ANGLES_FROM_TLT
            magnification=10000,
            samplingRate=self.samplingRate,
            voltage=300,
            dosePerFrame=2.1,
            tiltAxisAngle=87.2
        )
        protImportTS.setObjLabel('import tilt series')
        protImportTS = self.launchProtocol(protImportTS)
        tsSet = getattr(protImportTS, 'outputTiltSeries', None)

        # Validate output tomograms
        self.assertSetSize(tsSet, size=2)
        self.assertEqual(tsSet.getSamplingRate(), self.samplingRate)
        self.assertEqual(tsSet.getDim(), (1854, 1920, 1))

        return protImportTS

    def _tSCtfEstimateImod(self, protImportTS):
        print(magentaStr("\n==> Calculating the tilt series ctf:"))
        protTSCtfImod = self.newProtocol(
            ProtImodAutomaticCtfEstimation,
            inputSet=getattr(protImportTS, 'outputTiltSeries', None),
            angleRange=10
        )
        protTSCtfImod.setObjLabel('calculate TS CTF')
        protTSCtfImod = self.launchProtocol(protTSCtfImod)
        ctfSeriesSet = getattr(protTSCtfImod, 'outputSetOfCTFTomoSeries', None)

        # Validate output tomograms
        self.assertSetSize(ctfSeriesSet, size=2)

        return protTSCtfImod

    def _estimateCTF3D(self, protImportCoords3D, protTSCtfImod):
        print(magentaStr("\n==> Estimating the 3D CTF:"))
        protEstimateCTF3D = self.newProtocol(
            ProtRelionEstimateCTF3D,
            inputCoordinates=getattr(protImportCoords3D, 'outputCoordinates', None),
            inputSetCTFTomoSeries=getattr(protTSCtfImod, 'outputSetOfCTFTomoSeries', None),
            doseFilesPath=self.ds.getPath(),
            filesPattern='*ExpDose.txt',
            boxSize=72,
            ctf3dMode=CTF3D_PER_SUBVOLUME,
        )
        protEstimateCTF3D.setObjLabel('Estimate CTF 3D')
        protEstimateCTF3D = self.launchProtocol(protEstimateCTF3D)
        coord3DSet = getattr(protEstimateCTF3D, 'outputCoordinates', None)

        # Validate output tomograms
        self.assertSetSize(coord3DSet, size=150)
        self.assertEqual(coord3DSet.getSamplingRate(), self.samplingRate)
        self.assertEqual(coord3DSet.getBoxSize(), self.boxSize)
        # Output coordinates must have an attribute named _3dcftMrcFile, which stores the
        # path of the each ctf3D file
        for coord3d in coord3DSet:
            self.assertTrue(exists(coord3d._3dcftMrcFile.get()))

        return protEstimateCTF3D

    def _importVolume(self):
        print(magentaStr("\n==> Importing reference volume:"))
        protImportVol = self.newProtocol(
            ProtImportVolumes,
            emdbId=0,  # IMPORT_FROM_FILES
            filesPath=self.ds.getFile('refVol'),
            samplingRate=self.samplingRate,
        )
        protImportVol.setObjLabel('Import reference volume')
        protImportVol = self.launchProtocol(protImportVol)
        outVol = getattr(protImportVol, 'outputVolume', None)

        # Validate output tomograms
        self.assertEqual(outVol.getDim(), (self.boxSize, self.boxSize, self.boxSize))
        self.assertEqual(outVol.getSamplingRate(), self.samplingRate)

        return protImportVol

    def _importMask(self):
        print(magentaStr("\n==> Importing the mask:"))
        protImportMask = self.newProtocol(
            ProtImportMask,
            maskPath=self.ds.getFile('mask'),
            samplingRate=self.samplingRate,
        )
        protImportMask.setObjLabel('Import mask')
        protImportMask = self.launchProtocol(protImportMask)
        outVol = getattr(protImportMask, 'outputMask', None)

        # Validate output tomograms
        self.assertEqual(outVol.getDim(), (self.boxSize, self.boxSize, self.boxSize))
        self.assertEqual(outVol.getSamplingRate(), self.samplingRate)

        return protImportMask

    def _extractSubtomograms(self, protImportTomo, protEstimateCTF3D):
        print(magentaStr("\n==> Extracting subtomograms:"))
        protDynamoExtract = self.newProtocol(
            DynamoExtraction,
            inputCoordinates=getattr(protEstimateCTF3D, 'outputCoordinates', None),
            tomoSource=1,  # OTHER
            boxSize=self.boxSize,
            doInvert=True,
            inputTomograms=getattr(protImportTomo, 'outputTomograms', None)
        )
        protDynamoExtract.setObjLabel('extract subtomograms')
        protDynamoExtract = self.launchProtocol(protDynamoExtract)
        subtomoSet = getattr(protDynamoExtract, 'outputSetOfSubtomogram', None)

        # Validate output tomograms
        self.assertSetSize(subtomoSet, size=150)
        self.assertEqual(subtomoSet.getSamplingRate(), self.samplingRate)
        self.assertEqual(subtomoSet.getDim(), (self.boxSize, self.boxSize, self.boxSize))

        return protDynamoExtract

    def _classifySubtomograms(self, protExtractSubtomo, protImportRefVol, protImportMask):
        print(magentaStr("\n==> Classifying subtomograms:"))
        protClassifSubtomo = self.newProtocol(
            ProtRelionSubtomoClassif3D,
            threads=1,
            mpi=4,
            inputSubtomograms=getattr(protExtractSubtomo, 'outputSetOfSubtomogram', None),
            numberOfClasses=1,
            referenceVolume=getattr(protImportRefVol, 'outputVolume', None),
            doCTF=False,
            hasReferenceCTFCorrected=True,
            ctfMultiplied=True,
            referenceMask=getattr(protImportMask, 'outputMask', None),
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

        return protClassifSubtomo

    def _refineSubtomograms(self, protExtractSubtomo, protImportRefVol, protImportMask):
        print(magentaStr("\n==> Refining subtomograms:"))
        protRefineSubtomo = self.newProtocol(
            ProtRelionSubtomoRefine3D,
            threads=3,
            mpi=5,
            inputSubtomograms=getattr(protExtractSubtomo, 'outputSetOfSubtomogram', None),
            numberOfClasses=1,
            referenceVolume=getattr(protImportRefVol, 'outputVolume', None),
            doCTF=False,
            hasReferenceCTFCorrected=True,
            ctfMultiplied=True,
            referenceMask=getattr(protImportMask, 'outputMask', None),
            pooledSubtomos=10,
            doGpu=True,
            gpusToUse='0',
            numberOfIterations=5,
            extraParams='--sigma_tilt 3.667 --sigma_psi 3.667'
        )
        protRefineSubtomo.setObjLabel('refine subtomogram')
        protRefineSubtomo = self.launchProtocol(protRefineSubtomo)
        subTomoSet = getattr(protRefineSubtomo, 'outputParticles', None)
        avgSubTomo = getattr(protRefineSubtomo, 'outputVolume', None)

        # Validate output subtomograms
        self.assertSetSize(subTomoSet, size=150)
        self.assertEqual(subTomoSet.getSamplingRate(), self.samplingRate)
        self.assertEqual(subTomoSet.getDim(), (self.boxSize, self.boxSize, self.boxSize))

        # Validate output average subtomogram
        self.assertEqual(avgSubTomo.getSamplingRate(), self.samplingRate)
        self.assertEqual(avgSubTomo.getDim(), (self.boxSize, self.boxSize, self.boxSize))

        return protRefineSubtomo

    def _reconstructSubtomograms(self, protExtractSubtomo):
        print(magentaStr("\n==> Reconstructing subtomograms:"))
        protReconstructSubtomo = self.newProtocol(
            ProtRelionSubTomoReconstruct,
            threads=3,
            mpi=5,
            inputSubtomos=getattr(protExtractSubtomo, 'outputSetOfSubtomogram', None),
            maxRes=40
        )
        protReconstructSubtomo.setObjLabel('reconstruct subtomograms')
        protReconstructSubtomo = self.launchProtocol(protReconstructSubtomo)
        recTomo = getattr(protReconstructSubtomo, 'outputAvgSubtomogram', None)

        # Validate output average subtomograms
        self.assertEqual(recTomo.getSamplingRate(), self.samplingRate)
        self.assertEqual(recTomo.getDim(), (self.boxSize, self.boxSize, self.boxSize))

    def test_workflow(self):
        protImportTomo = self._importTomograms()
        protImportCoords3D = self._importCoordinates3D(protImportTomo)
        protImportTS = self._importTiltSeries()
        protTSCtfEstimationImod = self._tSCtfEstimateImod(protImportTS)
        protEstimateCTF3D = self._estimateCTF3D(protImportCoords3D, protTSCtfEstimationImod)
        protExtractSubtomo = self._extractSubtomograms(protImportTomo, protEstimateCTF3D)
        protImportRefVol = self._importVolume()
        protImportMask = self._importMask()
        self._classifySubtomograms(protExtractSubtomo, protImportRefVol, protImportMask)
        self._refineSubtomograms(protExtractSubtomo, protImportRefVol, protImportMask)
        self. _reconstructSubtomograms(protExtractSubtomo)
