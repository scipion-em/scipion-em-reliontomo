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
from pyworkflow.tests import BaseTest, setupTestProject, DataSet
from pyworkflow.utils import magentaStr
from reliontomo.protocols import ProtImportCoordinates3DFromStar, ProtImportSubtomogramsFromStar


class TestImportFromStarFile(BaseTest):

    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.dataset = DataSet.getDataSet('emd_10439')
        cls.samplingRate = 13.68
        cls.boxSize = 44
        cls.tomoId = 'emd_10439'

    def _runImportCoords3dFromStarFile(self, starFile, samplingRate=None):
        protImportCoords3dFromStar = self.newProtocol(ProtImportCoordinates3DFromStar,
                                                      starFile=starFile,
                                                      samplingRate=samplingRate,
                                                      boxSize=self.boxSize)

        self.launchProtocol(protImportCoords3dFromStar)
        self._checkCoordinatesAndTomograms(protImportCoords3dFromStar, coordSetSize=2339)

    def _checkCoordinatesAndTomograms(self, prot, coordSetSize=None):
        outputTomoSet = getattr(prot, 'outputTomograms', None)
        outputCoordsSet = getattr(prot, 'outputCoordinates', None)
        # Check the output set of 3D coordinates
        self.assertTrue(outputCoordsSet, 'No 3D coordinates were registered in the protocol output.')
        self.assertSetSize(outputCoordsSet, size=coordSetSize)
        self.assertEqual(outputCoordsSet.getBoxSize(), self.boxSize)
        self.assertEqual(outputCoordsSet.getSamplingRate(), self.samplingRate)
        self.assertTrue(outputCoordsSet.getPrecedents(), outputTomoSet)
        [self.assertEqual(coord.getTomoId(), self.tomoId) for coord in outputCoordsSet]

        # Check the output set of tomograms
        self.assertTrue(outputTomoSet, 'No tomograms were registered in the protocol output.')
        self.assertSetSize(outputTomoSet, size=1)
        self.assertEqual(outputTomoSet.getSamplingRate(), self.samplingRate)
        [self.assertEqual(tomo.getTsId(), self.tomoId) for tomo in outputTomoSet]

    def testImport3dCoordsFromStarFile_01(self):
        print(magentaStr("\n==> Importing coordinates 3D from a star file. Sampling rate read from protocol form:"))
        self._runImportCoords3dFromStarFile(self.dataset.getFile('coords3dStarFile'), self.samplingRate)

    def testImport3dCoordsFromStarFile_02(self):
        print(magentaStr("\n==> Importing coordinates 3D from a star file. Sampling rate read from the star file "
                         "(field rlnDetectorPixelSize):"))
        self._runImportCoords3dFromStarFile(self.dataset.getFile('coords3dStarFileWithSRate'))

    def _runImportSubtomogramsFromStarFile(self, starFile, samplingRate=None):
        protImportSubtomogramsFromStar = self.newProtocol(ProtImportSubtomogramsFromStar,
                                                          starFile=starFile,
                                                          samplingRate=samplingRate,
                                                          boxSize=self.boxSize)

        self.launchProtocol(protImportSubtomogramsFromStar)
        self._checkCoordinatesAndTomograms(protImportSubtomogramsFromStar, coordSetSize=7)
        self._checkSubtomograms(protImportSubtomogramsFromStar)

    def _checkSubtomograms(self, prot):
        subtomoSet = getattr(prot, 'outputSubtomograms', None)
        self.assertSetSize(subtomoSet, size=7)
        self.assertEqual(subtomoSet.getSamplingRate(), self.samplingRate)
        self.assertEqual(subtomoSet.getDim(), (self.boxSize, self.boxSize, self.boxSize))

    def testImportSubtomogramsFromStarFile(self):
        print(magentaStr("\n==> Importing subtomograms from a star file:"))
        self._runImportSubtomogramsFromStarFile(self.dataset.getFile('subtomogramsStarFile'), self.samplingRate)




