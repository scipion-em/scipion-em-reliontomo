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
from imod.protocols.protocol_base import OUTPUT_TOMOGRAMS_NAME
from pyworkflow.tests import BaseTest, setupTestProject, DataSet
from pyworkflow.utils import magentaStr, yellowStr
from reliontomo import Plugin
from reliontomo.protocols import ProtImportCoordinates3DFromStar
from reliontomo.protocols.protocol_base_import_from_star import importCoordsOutputs, IS_RE5_PICKING_ATTR
from reliontomo.protocols.protocol_import_subtomograms_from_star import outputObjects as importSubtomosOutputs
from tomo.constants import BOTTOM_LEFT_CORNER
from tomo.protocols import ProtImportTomograms
from tomo.protocols.protocol_import_tomograms import OUTPUT_NAME
from tomo.tests import EMD_10439, DataSetEmd10439, DataSetRe4STATuto, RE4_STA_TUTO, RE5_STA, DataSetRe5STA
from tomo.tests.test_base_centralized_layer import TestBaseCentralizedLayer

IS_RE_40 = Plugin.isRe40()
IS_RE_50 = Plugin.isRe50()


class TestImportFromStarFile(BaseTest):
    inTomoSet = None
    inTomoSetBin2 = None
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
        cls.inTomoSetBin2 = cls._normalizeTomograms()
        cls.coords1 = cls._runImportCoords3dFromStarFile(cls.dataset.getFile(DataSetEmd10439.coords3dStarFile.name),
                                                         samplingRate=cls.samplingRate,
                                                         inTomos=cls.inTomoSet)
        cls.coords2 = cls._runImportCoords3dFromStarFile(
            cls.dataset.getFile(DataSetEmd10439.coords3dStarFileWithSRate.name),
            inTomos=cls.inTomoSet)

    @classmethod
    def _importTomograms(cls):
        print(magentaStr("\n==> Importing data - tomograms:"))
        protImportTomogram = cls.newProtocol(ProtImportTomograms,
                                             filesPath=cls.dataset.getFile(DataSetEmd10439.tomoEmd10439.name),
                                             samplingRate=cls.samplingRate)

        cls.launchProtocol(protImportTomogram)
        outputTomos = protImportTomogram.Tomograms
        cls.assertIsNotNone(outputTomos, 'No tomograms were generated.')
        return outputTomos

    @classmethod
    def _normalizeTomograms(cls):
        # Generate a tomogram with a different sampling rate
        protNormTomogram = cls.newProtocol(ProtImodTomoNormalization,
                                           inputSetOfTomograms=cls.inTomoSet,
                                           binning=2)
        cls.launchProtocol(protNormTomogram)
        outputTomos = protNormTomogram.Tomograms
        cls.assertIsNotNone(outputTomos, 'No tomograms were generated.')
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
        cls.assertIsNotNone(outCoords, 'No coordinates were generated.')
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

        # Import the coordinates from a star file with the binned tomogram
        coords3 = self._runImportCoords3dFromStarFile(self.dataset.getFile(DataSetEmd10439.coords3dStarFile.name),
                                                      samplingRate=self.samplingRate,
                                                      inTomos=self.inTomoSetBin2)

        # These coordinates are referred to a tomogram which was binned a factor of 2 respecting the original tomogram
        # from this dataset, so the coordinates from one to another import protocols should present a scale factor of 2
        binFactor = 2
        self.assertEqual(coords3.getSamplingRate(), binFactor * self.coords1.getSamplingRate())
        for i in range(1, self.coordSetSize + 1):
            self.assertEqual(binFactor * coords3[i].getX(BOTTOM_LEFT_CORNER), self.coords1[i].getX(BOTTOM_LEFT_CORNER))
            self.assertEqual(binFactor * coords3[i].getY(BOTTOM_LEFT_CORNER), self.coords1[i].getY(BOTTOM_LEFT_CORNER))
            self.assertEqual(binFactor * coords3[i].getZ(BOTTOM_LEFT_CORNER), self.coords1[i].getZ(BOTTOM_LEFT_CORNER))

    def _runImportSubtomogramsFromStarFile(self, starFile, inTomoSet, binning=1):
        from reliontomo.protocols import ProtImportSubtomogramsFromStar
        protImportSubtomogramsFromStar = self.newProtocol(ProtImportSubtomogramsFromStar,
                                                          starFile=starFile,
                                                          inTomos=inTomoSet,
                                                          samplingRate=self.samplingRate,
                                                          boxSize=self.boxSize)

        self.launchProtocol(protImportSubtomogramsFromStar)
        outCoords = getattr(protImportSubtomogramsFromStar, importSubtomosOutputs.coordinates.name)
        outSubtomos = getattr(protImportSubtomogramsFromStar, importSubtomosOutputs.subtomograms.name)
        self._checkCoordinates(outCoords, coordSetSize=7, inTomos=inTomoSet)
        self._checkSubtomograms(outSubtomos, binning=binning)

    def _checkSubtomograms(self, subtomoSet, binning):
        self.assertSetSize(subtomoSet, size=7)
        self.assertEqual(subtomoSet.getSamplingRate(), self.samplingRate * binning)
        self.assertEqual(subtomoSet.getDim(), (self.boxSize, self.boxSize, self.boxSize))

    def testImportSubtomogramsFromStarFile_01(self):
        if IS_RE_40:
            print(magentaStr("\n==> Importing subtomograms from a star file:"))
            self._runImportSubtomogramsFromStarFile(self.dataset.getFile(DataSetEmd10439.subtomogramsStarFile.name),
                                                    self.inTomoSet)
        else:
            print(yellowStr('Relion 5 detected. Test for protocol "Import subtomograms from star file" skipped.'))

    def testImportSubtomogramsFromStarFile_02(self):
        if IS_RE_40:
            print(magentaStr("\n==> Importing subtomograms from a star file with a sampling rate different than the "
                             "tomograms sampling rate:"))
            self._runImportSubtomogramsFromStarFile(self.dataset.getFile(DataSetEmd10439.subtomogramsStarFile.name),
                                                    self.inTomoSetBin2,
                                                    binning=2)
        else:
            print(yellowStr('Relion 5 detected. Test for protocol "Import subtomograms from star file" skipped.'))


class TestRelion5ImportFromStarFile(TestBaseCentralizedLayer):
    bin4 = 4
    boxSizeBin4 = DataSetRe4STATuto.boxSizeBin4.value

    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.dsRe4 = DataSet.getDataSet(RE4_STA_TUTO)
        cls.tomogramsImported = cls._importTomograms()

    @classmethod
    def _importTomograms(cls):
        print(magentaStr("\n==> Importing data - tomograms:"))
        protImportTomogram = cls.newProtocol(ProtImportTomograms,
                                             filesPath=cls.dsRe4.getFile(DataSetRe4STATuto.tomosPath.name),
                                             filesPattern=DataSetRe4STATuto.tomosPattern.value,
                                             samplingRate=DataSetRe4STATuto.sRateBin4.value)

        cls.launchProtocol(protImportTomogram)
        outputTomos = getattr(protImportTomogram, OUTPUT_NAME, None)
        cls.assertIsNotNone(outputTomos, 'No tomograms were generated.')
        return outputTomos

    @classmethod
    def _binTomograms(cls):
        print(magentaStr("\n==> Binning the tomograms:"))
        protNormTomogram = cls.newProtocol(ProtImodTomoNormalization,
                                           inputSetOfTomograms=cls.tomogramsImported,
                                           binning=2)
        cls.launchProtocol(protNormTomogram)
        outputTomos = getattr(protNormTomogram, OUTPUT_TOMOGRAMS_NAME, None)
        cls.assertIsNotNone(outputTomos, 'No tomograms were generated.')
        return outputTomos

    @classmethod
    def _runImportCoords3dFromStar(cls, starFile, samplingRate=None, inTomos=None, boxSize=None):
        protImportCoords3dFromStar = cls.newProtocol(ProtImportCoordinates3DFromStar,
                                                     starFile=starFile,
                                                     inTomos=inTomos,
                                                     samplingRate=samplingRate,
                                                     boxSize=boxSize)

        cls.launchProtocol(protImportCoords3dFromStar)
        outCoords = getattr(protImportCoords3dFromStar, importCoordsOutputs.coordinates.name)
        cls.assertIsNotNone(outCoords, 'No coordinates were generated.')
        return outCoords

    def _runTest(self, tomograms=None, sRate=None, boxSize=None, isRelion5Picking=None):
        starFile = self.dsRe4.getFile(DataSetRe4STATuto.coordsStarSubset.name)
        importedCoords = self._runImportCoords3dFromStar(starFile,
                                                         inTomos=tomograms,
                                                         samplingRate=sRate,
                                                         boxSize=boxSize)
        # Check the results
        self.checkCoordinates(importedCoords,
                              expectedSetSize=DataSetRe4STATuto.nCoordsTotal.value,
                              expectedBoxSize=boxSize,
                              expectedSRate=tomograms.getSamplingRate(),  # The coords are scaled to the
                              # size of the tomograms introduced
                              orientedParticles=True)  # Oriented picking
        self.assertEqual(getattr(importedCoords, IS_RE5_PICKING_ATTR, -1), isRelion5Picking)

    def testImportCoords_01(self):
        print(magentaStr("\n==> Importing coordinates 3D from a star file."
                         "\n\t- The coordinates sampling rate is read from protocol form."))
        tomograms = self.tomogramsImported
        sRate = DataSetRe4STATuto.unbinnedPixSize.value  # Coords provided by teh tutorial of Relion 4 are at bin 1
        boxSize = self.boxSizeBin4
        self._runTest(tomograms=tomograms, sRate=sRate, boxSize=boxSize, isRelion5Picking=False)

    def testImportCoords_02(self):
        print(magentaStr("\n==> Importing coordinates 3D from a star file."
                         "\n\t- The coordinates sampling rate sampling rate is not provided, so the one from the "
                         "tomograms is assumed."))
        tomograms = self.tomogramsImported
        boxSize = self.boxSizeBin4
        self._runTest(tomograms=tomograms, boxSize=boxSize, isRelion5Picking=False)

    def testImportCoords_03(self):
        print(magentaStr("\n==> Importing coordinates 3D from a star file."
                         "\n\t- The tomograms introduced are binned a factor of 2."))
        tomograms = self._binTomograms()
        sRate = DataSetRe4STATuto.unbinnedPixSize.value  # Coords provided by teh tutorial of Relion 4 are at bin 1
        boxSize = self.boxSizeBin4 / 2
        self._runTest(tomograms=tomograms, sRate=sRate, boxSize=boxSize, isRelion5Picking=False)


class TestImportRe5NativeCoordsFromStarFile(TestBaseCentralizedLayer):
    """Import coordinates picked with native Relion 5"""

    @classmethod
    def setUpClass(cls):
        if IS_RE_50:
            setupTestProject(cls)
            cls.ds = DataSet.getDataSet(RE5_STA)
            cls.importedTomos = cls._importTomograms()

    @classmethod
    def _importTomograms(cls):
        print(magentaStr("\n==> Importing data - tomograms:"))
        protImportTomogram = cls.newProtocol(ProtImportTomograms,
                                             filesPath=cls.ds.getFile(DataSetRe5STA.tomosDir.name),
                                             filesPattern='*',
                                             samplingRate=DataSetRe4STATuto.sRateBin4.value)

        cls.launchProtocol(protImportTomogram)
        outputTomos = getattr(protImportTomogram, OUTPUT_NAME, None)
        cls.assertIsNotNone(outputTomos, 'No tomograms were generated.')
        return outputTomos

    @classmethod
    def _runImportCoords3dFromStar(cls, starFile, samplingRate=None, inTomos=None, boxSize=None):
        protImportCoords3dFromStar = cls.newProtocol(ProtImportCoordinates3DFromStar,
                                                     starFile=starFile,
                                                     inTomos=inTomos,
                                                     samplingRate=samplingRate,
                                                     boxSize=boxSize)

        cls.launchProtocol(protImportCoords3dFromStar)
        outCoords = getattr(protImportCoords3dFromStar, importCoordsOutputs.coordinates.name)
        cls.assertIsNotNone(outCoords, 'No coordinates were generated.')
        return outCoords

    def _runTest(self, coordsSRate=None):
        if IS_RE_50:
            sRateMsg = '' if coordsSRate else 'not'
            print(magentaStr(f"\n==> Import coordinates picked with Relion 5:"
                             f"\n\t- Coordinates sampling rate {sRateMsg} provided."))
            starFile = self.ds.getFile(DataSetRe5STA.coordsPickedWithRe5Star.name)
            boxSize = DataSetRe5STA.boxSize.value
            tomograms = self.importedTomos
            importedCoords = self._runImportCoords3dFromStar(starFile,
                                                             inTomos=tomograms,
                                                             samplingRate=coordsSRate,
                                                             boxSize=boxSize)
            # Check the results
            isRelion5Picking = True
            self.checkCoordinates(importedCoords,
                                  expectedSetSize=DataSetRe5STA.nCoords.value,
                                  expectedBoxSize=boxSize,
                                  expectedSRate=tomograms.getSamplingRate(),  # The coords are scaled to the
                                  # size of the tomograms introduced
                                  orientedParticles=True)  # Oriented picking
            self.assertEqual(getattr(importedCoords, IS_RE5_PICKING_ATTR, -1), isRelion5Picking)
        else:
            print(yellowStr('Relion 4 detected. Test for protocol "Import coordinates picked with Relion 5" skipped.'))

    def testImportsCoordsNativeRe5_01(self):
        self._runTest(coordsSRate=DataSetRe5STA.coordsSRate.value)

    def testImportsCoordsNativeRe5_02(self):
        self._runTest()
