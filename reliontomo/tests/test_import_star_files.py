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
import reliontomo
from pyworkflow.tests import BaseTest, setupTestProject, DataSet
from pyworkflow.utils import magentaStr


class TestTomoImportSetOfCoordinates3D(BaseTest):
    """This class check if the protocol to import set of coordinates 3d works properly."""

    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.dataset = DataSet.getDataSet('reliontomo')
        cls.samplingRate = 13.68

    def _runImportCoords3dFromStarFile(self, starFile, samplingRate=None):
        boxSize = 44
        tomoId = 'emd_10439'
        protImportCoords3dFromStar = self.newProtocol(reliontomo.protocols.ProtImportCoordinates3DFromStar,
                                                      starFile=starFile,
                                                      samplingRate=samplingRate,
                                                      boxSize=boxSize)

        self.launchProtocol(protImportCoords3dFromStar)
        outputTomoSet = getattr(protImportCoords3dFromStar, 'outputSetOfTomograms', None)
        outputCoordsSet = getattr(protImportCoords3dFromStar, 'outputCoordinates', None)

        # Check the output set of 3D coordinates
        self.assertTrue(outputCoordsSet, 'No 3D coordinates were registered in the protocol output.')
        self.assertSetSize(outputCoordsSet, size=2339)
        self.assertEqual(outputCoordsSet.getBoxSize(), boxSize)
        self.assertEqual(outputCoordsSet.getSamplingRate(), self.samplingRate)
        self.assertTrue(outputCoordsSet.getPrecedents(), outputTomoSet)
        [self.assertEqual(coord.getTomoId(), tomoId) for coord in outputCoordsSet]

        # Check the output set of tomograms
        self.assertTrue(outputTomoSet, 'No tomograms were registered in the protocol output.')
        self.assertSetSize(outputTomoSet, size=1)
        self.assertEqual(outputTomoSet.getSamplingRate(), self.samplingRate)
        [self.assertEqual(tomo.getTsId(), tomoId) for tomo in outputTomoSet]

    def testImport3dCoordsFromStarFile_01(self):
        print(magentaStr("\n==> Importing coordinates 3D from a star file. Sampling rate read from protocol form:"))
        self._runImportCoords3dFromStarFile(self.dataset.getFile('coords3dStarFile'), self.samplingRate)

    def testImport3dCoordsFromStarFile_02(self):
        print(magentaStr("\n==> Importing coordinates 3D from a star file. Sampling rate read from the star file "
                         "(field rlnDetectorPixelSize):"))
        self._runImportCoords3dFromStarFile(self.dataset.getFile('coords3dStarFileWithSRate'))

