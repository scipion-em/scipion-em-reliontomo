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
from pyworkflow.tests import BaseTest, DataSet
from pyworkflow.utils import magentaStr
from reliontomo.constants import STAR_FILES_EQUAL, STAR_DIFF_SIZE, STAR_DIFF_LABELS, STAR_DIFF_VALUES, OPTICS_TABLE, \
    PARTICLES_TABLE, GLOBAL_TABLE, TOMO_PARTICLE_ID, SUBTOMO_NAME, CTF_IMAGE, CLASS_NUMBER, RANDOM_SUBSET, \
    TILT_SERIES_NAME
from reliontomo.objects import StarFileComparer
from reliontomo.tests import RE4_TOMO, DataSetRe4Tomo

# Tomograms star tables
TS_45_TABLE = 'TS_45'
TS_54_TABLE = 'TS_54'


class TestStarFileComparer(BaseTest):

    @classmethod
    def setUpClass(cls):
        cls.dataset = DataSet.getDataSet(RE4_TOMO)

    def test_01_checkOKStar(self):
        """Compare a star file with itself."""
        starFile1 = self.dataset.getFile(DataSetRe4Tomo.okStarFile.name)
        filesEqual = True
        diffSize = False
        diffLabels = False
        diffValues = False
        testName = self.test_01_checkOKStar.__name__
        # Check optics table
        sfc = StarFileComparer(starFile1, starFile1, OPTICS_TABLE)
        self._checkComparison(sfc.compare(),
                              testName=testName,
                              filesEqual=filesEqual,
                              diffSize=diffSize,
                              diffLabels=diffLabels,
                              diffValues=diffValues)
        # Check particles table
        sfc = StarFileComparer(starFile1, starFile1, PARTICLES_TABLE)
        self._checkComparison(sfc.compare(),
                              testName=testName,
                              filesEqual=filesEqual,
                              diffSize=diffSize,
                              diffLabels=diffLabels,
                              diffValues=diffValues)

    def test_02_checkDiffNumberOfRows(self):
        """Compare an OK star file with a copy of itself excepting 5 lines removed in table named particles."""
        starFile1 = self.dataset.getFile(DataSetRe4Tomo.okStarFile.name)
        starFile2 = self.dataset.getFile(DataSetRe4Tomo.diffSizeStar.name)
        testName = self.test_02_checkDiffNumberOfRows.__name__
        # Check optics table
        sfc = StarFileComparer(starFile1, starFile2, OPTICS_TABLE)
        self._checkComparison(sfc.compare(),
                              testName=testName,
                              filesEqual=True,
                              diffSize=False,
                              diffLabels=False,
                              diffValues=False)
        # Check particles table
        sfc = StarFileComparer(starFile1, starFile2, PARTICLES_TABLE)
        compareMsg = sfc.compare()
        self._checkComparison(compareMsg,
                              testName=testName,
                              filesEqual=False,
                              diffSize=True,
                              diffLabels=False,
                              diffValues=False)
        self.assertTrue('108 != 103' in compareMsg)

    def test_03_checkDiffLabels(self):
        """Compare an OK star file with a copy of itself excepting 2 labels removed in table named particles,
        concretely the labels rlnClassNumber and rlnRandomSubset."""
        starFile1 = self.dataset.getFile(DataSetRe4Tomo.okStarFile.name)
        starFile2 = self.dataset.getFile(DataSetRe4Tomo.diffLabelsStar.name)
        testName = self.test_03_checkDiffLabels.__name__
        # Check optics table
        sfc = StarFileComparer(starFile1, starFile2, OPTICS_TABLE)
        self._checkComparison(sfc.compare(),
                              testName=testName,
                              filesEqual=True,
                              diffSize=False,
                              diffLabels=False,
                              diffValues=False)
        # Check particles table
        sfc = StarFileComparer(starFile1, starFile2, PARTICLES_TABLE)
        strings2check = ['DIFF LABELS', CLASS_NUMBER, RANDOM_SUBSET]
        compareMsg = sfc.compare()
        self._checkComparison(compareMsg,
                              testName=testName,
                              filesEqual=False,
                              diffSize=False,
                              diffLabels=True,
                              diffValues=False,
                              strings2check=strings2check)

    def test_04_checkDiffValues(self):
        """Compare an OK star file with a copy of itself excepting some values changed."""
        starFile1 = self.dataset.getFile(DataSetRe4Tomo.okStarFile.name)
        starFile2 = self.dataset.getFile(DataSetRe4Tomo.diffValuesStar.name)
        testName = self.test_04_checkDiffValues.__name__
        # Check optics table
        sfc = StarFileComparer(starFile1, starFile2, OPTICS_TABLE)
        self._checkComparison(sfc.compare(),
                              testName=testName,
                              filesEqual=True,
                              diffSize=False,
                              diffLabels=False,
                              diffValues=False)
        # Check particles table
        strings2check = ['ROW 2', 'ROW 12', 'ROW 41', 'rlnCoordinateX', 'rlnCoordinateX', '664.0 != 665.0']
        sfc = StarFileComparer(starFile1, starFile2, PARTICLES_TABLE)
        self._checkComparison(sfc.compare(),
                              testName=testName,
                              filesEqual=False,
                              diffSize=False,
                              diffLabels=False,
                              diffValues=True,
                              strings2check=strings2check)

    def test_05_prepareTomosInVsOut(self):
        starFile1 = self.dataset.getFile(DataSetRe4Tomo.prepareTomosStarScipion.name)
        starFile2 = self.dataset.getFile(DataSetRe4Tomo.prepareTomosStarRelion.name)
        tablesList = [GLOBAL_TABLE, TS_45_TABLE, TS_54_TABLE]
        testName = self.test_05_prepareTomosInVsOut.__name__
        excludeLabelsList = [TILT_SERIES_NAME]  # paths may differ
        self._checkMultipleTablesOk(starFile1, starFile2, tablesList,
                                    testName=testName,
                                    excludeLabelsList=excludeLabelsList)

    def test_06_prepareParticlesInVsOut(self):
        starFile1 = self.dataset.getFile(DataSetRe4Tomo.preparePartcilesStarScipion.name)
        starFile2 = self.dataset.getFile(DataSetRe4Tomo.prepareParticlesStarRelion.name)
        tablesList = [OPTICS_TABLE, PARTICLES_TABLE]
        testName = self.test_06_prepareParticlesInVsOut.__name__
        excludeLabelsList = ['rlnTomoParticleId']
        self._checkMultipleTablesOk(starFile1, starFile2, tablesList,
                                    testName=testName,
                                    excludeLabelsList=excludeLabelsList)

    def test_07_makePSubtomosInVsOut(self):
        starFile1 = self.dataset.getFile(DataSetRe4Tomo.makePSubtomosStarScipion.name)
        starFile2 = self.dataset.getFile(DataSetRe4Tomo.makePSubtomosStarRelion.name)
        tablesList = [OPTICS_TABLE, PARTICLES_TABLE]
        testName = self.test_07_makePSubtomosInVsOut.__name__
        excludeLabelsList = [TOMO_PARTICLE_ID, SUBTOMO_NAME, CTF_IMAGE]  # Paths may be different, but tsId are enough for the match
        self._checkMultipleTablesOk(starFile1, starFile2, tablesList,
                                    testName=testName,
                                    excludeLabelsList=excludeLabelsList)

    def test_08_makePSubtomosInVsOut(self):
        starFile1 = self.dataset.getFile(DataSetRe4Tomo.zShiftedScipion.name)
        starFile2 = self.dataset.getFile(DataSetRe4Tomo.zShiftedRelion.name)
        tablesList = [OPTICS_TABLE, PARTICLES_TABLE]
        testName = self.test_07_makePSubtomosInVsOut.__name__
        excludeLabelsList = ['rlnTomoParticleId']
        self._checkMultipleTablesOk(starFile1, starFile2, tablesList,
                                    testName=testName,
                                    excludeLabelsList=excludeLabelsList)

    def _checkComparison(self, compareMsg, testName=None, filesEqual=None, diffSize=None,
                         diffLabels=None, diffValues=None, strings2check=None):
        if testName:
            print(magentaStr("\n==> %s   %s" % (testName, compareMsg)))
        self.assertEqual(STAR_FILES_EQUAL in compareMsg, filesEqual)
        self.assertEqual(STAR_DIFF_SIZE in compareMsg, diffSize)
        self.assertEqual(STAR_DIFF_LABELS in compareMsg, diffLabels)
        self.assertEqual(STAR_DIFF_VALUES in compareMsg, diffValues)
        if strings2check:
            [self.assertTrue(testStr in compareMsg) for testStr in strings2check]

    def _checkMultipleTablesOk(self, starFile1, starFile2, tablesList, testName=None, excludeLabelsList=None):
        for tableName in tablesList:
            sfc = StarFileComparer(starFile1, starFile2, tableName)
            self._checkComparison(sfc.compare(excludeLabelsList=excludeLabelsList),
                                  testName=testName,
                                  filesEqual=True,
                                  diffSize=False,
                                  diffLabels=False,
                                  diffValues=False)
                                  # strings2check=strings2check)


