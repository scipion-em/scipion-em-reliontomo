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
from imod.protocols import ProtImodTomoNormalization
from pyworkflow.tests import BaseTest, setupTestProject, DataSet
from pyworkflow.utils import magentaStr
from reliontomo.protocols import ProtImportCoordinates3DFromStar, ProtImportSubtomogramsFromStar
from reliontomo.protocols.protocol_import_subtomograms_from_star import outputObjects as importSubtomosOutputs
from reliontomo.tests import DataSetEmd10439, EMD_10439, OUTPUT_TOMOS
from tomo.constants import BOTTOM_LEFT_CORNER
from tomo.protocols import ProtImportTomograms


class TestImportFromStarFile(BaseTest):

    inTomoSet = None
    boxSize = None
    samplingRate = None
    dataset = None

    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.dataset = DataSet.getDataSet(EMD_10439)
        cls.samplingRate = 13.68
        cls.boxSize = 44
        cls.coordSetSize = 2339
        cls.tomoId = EMD_10439
        cls.inTomoSet = cls._importTomograms()
        cls.coords1 = cls._runImportCoords3dFromStarFile(cls.dataset.getFile(DataSetEmd10439.coords3dStarFile.name),
                                                         samplingRate=cls.samplingRate,
                                                         inTomos=cls.inTomoSet)
        cls.coords2 = cls._runImportCoords3dFromStarFile(cls.dataset.getFile(DataSetEmd10439.coords3dStarFileWithSRate.name),
                                                         inTomos=cls.inTomoSet)

    @classmethod
    def _importTomograms(cls):
        print(magentaStr("\n==> Importing data - tomograms:"))
        protImportTomogram = cls.newProtocol(ProtImportTomograms,
                                             filesPath=cls.dataset.getFile(DataSetEmd10439.tomoEmd10439.name),
                                             samplingRate=cls.samplingRate)

        cls.launchProtocol(protImportTomogram)
        outputTomos = getattr(protImportTomogram, OUTPUT_TOMOS, None)
        cls.assertIsNotNone(outputTomos, 'No tomograms were genetated.')
        return outputTomos

    @classmethod
    def _runImportCoords3dFromStarFile(cls, starFile, samplingRate=None, inTomos=None):
        protImportCoords3dFromStar = cls.newProtocol(ProtImportCoordinates3DFromStar,
                                                     starFile=starFile,
                                                     inTomos=inTomos,
                                                     samplingRate=samplingRate,
                                                     boxSize=cls.boxSize)

        cls.launchProtocol(protImportCoords3dFromStar)
        outCoords = getattr(protImportCoords3dFromStar, importSubtomosOutputs.coordinates.name)
        cls.assertIsNotNone(outCoords, 'No coordinates were genetated.')
        return outCoords

    def _checkCoordinates(self, outputCoordsSet, coordSetSize=None, inTomos=None):
        # Check the output set of 3D coordinates
        self.assertTrue(outputCoordsSet, 'No 3D coordinates were registered in the protocol output.')
        self.assertSetSize(outputCoordsSet, size=coordSetSize)
        self.assertEqual(outputCoordsSet.getBoxSize(), self.boxSize)
        self.assertEqual(outputCoordsSet.getSamplingRate(), inTomos.getSamplingRate())
        self.assertTrue(outputCoordsSet.getPrecedents(), inTomos)
        [self.assertEqual(coord.getTomoId(), self.tomoId) for coord in outputCoordsSet]

    def testImport3dCoordsFromStarFile_01(self):
        print(magentaStr("\n==> Importing coordinates 3D from a star file. Sampling rate read from protocol form:"))
        self._checkCoordinates(self.coords1, coordSetSize=2339, inTomos=self.inTomoSet)

    def testImport3dCoordsFromStarFile_02(self):
        print(magentaStr("\n==> Importing coordinates 3D from a star file. Sampling rate read from the star file "
                         "(field rlnDetectorPixelSize):"))
        self._checkCoordinates(self.coords2, coordSetSize=2339, inTomos=self.inTomoSet)

    def testImport3dCoordsFromStarFile_03(self):
        print(magentaStr("\n==> Importing coordinates 3D from a star file with a sampling rate different than the "
                         "tomograms sampling rate:"))

        # Generate a tomogram with a different sampling rate
        protNormTomogram = self.newProtocol(ProtImodTomoNormalization,
                                            inputSetOfTomograms=self.inTomoSet,
                                            binning=2)
        self.launchProtocol(protNormTomogram)
        outputTomos = getattr(protNormTomogram, 'outputSetOfTomograms', None)
        self.assertIsNotNone(outputTomos, 'No tomograms were genetated.')

        # Import the coordinates from a star file with the binned tomogram
        coords3 = self._runImportCoords3dFromStarFile(self.dataset.getFile(DataSetEmd10439.coords3dStarFile.name),
                                                      samplingRate=self.samplingRate,
                                                      inTomos=outputTomos)

        # These coordinates are referred to a tomogram which was binned a factor of 2 respecting the original tomogram
        # from this dataset, so the coordinates from one to another import protocols should present a scale factor of 2
        binFactor = 2
        self.assertEqual(coords3.getSamplingRate(), binFactor * self.coords1.getSamplingRate())
        for i in range(1, self.coordSetSize + 1):
            self.assertEqual(binFactor * coords3[i].getX(BOTTOM_LEFT_CORNER), self.coords1[i].getX(BOTTOM_LEFT_CORNER))
            self.assertEqual(binFactor * coords3[i].getY(BOTTOM_LEFT_CORNER), self.coords1[i].getY(BOTTOM_LEFT_CORNER))
            self.assertEqual(binFactor * coords3[i].getZ(BOTTOM_LEFT_CORNER), self.coords1[i].getZ(BOTTOM_LEFT_CORNER))

    def _runImportSubtomogramsFromStarFile(self, starFile, samplingRate=None):
        protImportSubtomogramsFromStar = self.newProtocol(ProtImportSubtomogramsFromStar,
                                                          starFile=starFile,
                                                          inTomos=self.inTomoSet,
                                                          samplingRate=samplingRate,
                                                          boxSize=self.boxSize)

        self.launchProtocol(protImportSubtomogramsFromStar)
        outCoords = getattr(protImportSubtomogramsFromStar, importSubtomosOutputs.coordinates.name)
        outSubtomos = getattr(protImportSubtomogramsFromStar, importSubtomosOutputs.subtomograms.name)
        self._checkCoordinates(outCoords, coordSetSize=7, inTomos=self.inTomoSet)
        self._checkSubtomograms(outSubtomos)

    def _checkSubtomograms(self, subtomoSet):
        self.assertSetSize(subtomoSet, size=7)
        self.assertEqual(subtomoSet.getSamplingRate(), self.samplingRate)
        self.assertEqual(subtomoSet.getDim(), (self.boxSize, self.boxSize, self.boxSize))

    def testImportSubtomogramsFromStarFile(self):
        print(magentaStr("\n==> Importing subtomograms from a star file:"))
        self._runImportSubtomogramsFromStarFile(self.dataset.getFile(DataSetEmd10439.subtomogramsStarFile.name),
                                                self.inTomoSet.getSamplingRate())
