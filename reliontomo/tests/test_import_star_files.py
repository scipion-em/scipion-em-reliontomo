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


class TestTomoImportSetOfCoordinates3D(BaseTest):
    """This class check if the protocol to import set of coordinates 3d works properly."""

    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.dataset = DataSet.getDataSet('reliontomo')

    def testImportCoords3dFromStarFile(self):
        boxSize = 44
        samplingRate = 13.68
        protImportCoords3dFromStar = self.newProtocol(reliontomo.protocols.ProtImportCoordinates3DFromStar,
                                                      starFile=self.dataset.getFile('coords3dStarFile'),
                                                      samplingRate=samplingRate,
                                                      boxSize=boxSize)

        self.launchProtocol(protImportCoords3dFromStar)
        outputTomoSet = getattr(protImportCoords3dFromStar, 'outputSetOfTomograms', None)
        outputCoordsSet = getattr(protImportCoords3dFromStar, 'outputCoordinates', None)

        # Check output set of coordinates
        self.assertTrue(outputCoordsSet, 'There was a problem with coordinates 3d output')
        self.assertSetSize(outputCoordsSet, size=2339)
        self.assertEqual(outputCoordsSet.getBoxSize(), boxSize)
        self.assertEqual(outputCoordsSet.getSamplingRate(), samplingRate)
        self.assertTrue(outputCoordsSet.getPrecedents(), outputTomoSet)
