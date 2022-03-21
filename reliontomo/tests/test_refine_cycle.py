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
from reliontomo.protocols import ProtImportCoordinates3DFromStar
from reliontomo.tests import RE4_TOMO, DataSetRe4Tomo, OUTPUT_TOMOS, OUTPUT_COORDS
from tomo.protocols import ProtImportTomograms, ProtImportTs


class TestImportFromStarFile(BaseTest):
    inTomoSet = None
    boxSize = None
    dataset = None
    samplingRate = None

    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.samplingRate = 1.35
        cls.boxSize = 96
        cls.dataset = DataSet.getDataSet(RE4_TOMO)
        cls.inTomoSet = cls._importTomograms()
        cls.inCoords = cls._runImportCoords3dFromStarFile()

    @classmethod
    def _importTomograms(cls):
        print(magentaStr("\n==> Importing data - tomograms:"))
        tomogramsBinning = 4
        protImportTomogram = cls.newProtocol(ProtImportTomograms,
                                             filesPath=cls.dataset.getFile(DataSetRe4Tomo.tomogram.name),
                                             samplingRate=tomogramsBinning * cls.samplingRate)

        cls.launchProtocol(protImportTomogram)
        outputTomos = getattr(protImportTomogram, OUTPUT_TOMOS, None)
        cls.assertIsNotNone(outputTomos, 'No tomograms were genetated.')
        return outputTomos

    @classmethod
    def _runImportCoords3dFromStarFile(cls):
        protImportCoords3dFromStar = cls.newProtocol(ProtImportCoordinates3DFromStar,
                                                     starFile=cls.dataset.getFile(DataSetRe4Tomo.coordinates.name),
                                                     inTomos=cls.inTomoSet,
                                                     samplingRate=1.35,
                                                     boxSize=cls.boxSize)

        cls.launchProtocol(protImportCoords3dFromStar)
        outCoords = getattr(protImportCoords3dFromStar, OUTPUT_COORDS, None)
        cls.assertIsNotNone(outCoords, 'No coordinates were genetated.')
        return outCoords

    @classmethod
    def _runImportTS(cls):
        protImportTS = cls.newProtocol(ProtImportTs,
                                       filesPath=cls.dataset.getFile(DataSetRe4Tomo.eTomoDir),
                                       filesPattern='*/*.mdoc')
        cls.launchProtocol(protImportTS)
        outputTS = getattr(protImportTS, 'outputTiltSeries', None)
        cls.assertIsNotNone(outputTS, 'No tilt series were genetated.')
        return outputTS

    @classmethod
    def _estimateCTF(cls):
        pass
