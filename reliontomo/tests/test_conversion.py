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
from reliontomo.convert.convertBase import getTransformInfoFromCoordOrSubtomo, genTransformMatrix
from pyworkflow.tests import BaseTest
from pwem.objects import Transform
from pwem.convert.transformations import euler_matrix
import numpy as np
from reliontomo.objects import RelionPSubtomogram


class TestTransformationConversion(BaseTest):


    def test_transformations(self):
        """Test conversion o alignment information from and to relion"""
        self._test_set_transformations()

    def getInitialValues(self):
        """ Returns initial values, Use commented code to test a specific case"""
        fromRealCase ="3.290925     -1.08901    40.351380    65.410000   175.090000    18.310000"
        fromRealCase = fromRealCase.split()
        return float(fromRealCase[0]), float(fromRealCase[1]), float(fromRealCase[2]), \
            float(fromRealCase[3]), float(fromRealCase[4]), float(fromRealCase[5]),
    def _test_set_transformations(self):


        xShift ,yShift ,zShift, rot, tilt, psi =self.getInitialValues()

        self._test_single_conversion(psi, rot, tilt, 0, 0, 0, "ONLY ANGLES")
        self._test_single_conversion(0, 0, 0, xShift, yShift, zShift, "ONLY SHIFTS")
        self._test_single_conversion(psi, rot, tilt, xShift, yShift, zShift, "COMPLETE")

    def _test_single_conversion(self, psi, rot, tilt, xShift, yShift, zShift, testName):

        self._test_direct_conversion(psi, rot, tilt, xShift, yShift, zShift, testName)

        rotRad = np.deg2rad(rot)
        titlRad = np.deg2rad(tilt)
        psiRad = np.deg2rad(psi)
        alignmentMatrix = euler_matrix(rotRad, titlRad, psiRad, axes="szyz")
        alignmentMatrix[0, 3] = xShift
        alignmentMatrix[1, 3] = yShift
        alignmentMatrix[2, 3] = zShift
        # Convert to relion
        # We need a faked subtomogram or coordinate
        subtomo = RelionPSubtomogram()
        t = Transform(matrix=alignmentMatrix)
        subtomo.setTransform(t)
        rlnAngles, rlnShifts = getTransformInfoFromCoordOrSubtomo(subtomo, samplingRate=2)
        returningMatrix = genTransformMatrix(rlnShifts[0], rlnShifts[1], rlnShifts[2],
                                             rlnAngles[0], rlnAngles[1], rlnAngles[2], 2)

        ok = np.allclose(alignmentMatrix, returningMatrix)

        if not ok:
            print("*******", testName, "*******")
            print("EXPECTED MATRIX:\n %s" % np.array_str(alignmentMatrix, precision=2, suppress_small=True))
            print("RESULT MATRIX:\n %s" % np.array_str(returningMatrix, precision=2, suppress_small=True))

        self.assertTrue(ok, "Relion-Scipion transformations %s is wrong." % (testName))


    def _test_direct_conversion(self, psi, rot, tilt, xShift, yShift, zShift, testName):

        # Matrix for scipion
        alignmentMatrix = genTransformMatrix(xShift, yShift, zShift, rot, tilt, psi, 1.35)

        subtomo = RelionPSubtomogram()
        t = Transform(matrix=alignmentMatrix)
        subtomo.setTransform(t)

        # From scipion to relion
        rlnAngles, rlnShifts = getTransformInfoFromCoordOrSubtomo(subtomo, samplingRate=1.35)

        self.assertAlmostEqual(rlnAngles[0], rot, places=2,msg="Rot not expected for direct %s" % testName)
        self.assertAlmostEqual(rlnAngles[1], tilt, places=2,msg="Tilt not expected for direct %s" % testName)
        self.assertAlmostEqual(rlnAngles[2], psi, places=2,msg="Psi not expected for direct %s" % testName)

        self.assertAlmostEqual(rlnShifts[0], xShift, places=2,msg="xShift not expected for direct %s" % testName)
        self.assertAlmostEqual(rlnShifts[1], yShift, places=2,msg="yShift not expected for direct %s" % testName)
        self.assertAlmostEqual(rlnShifts[2], zShift, places=2,msg="zShift not expected for direct %s" % testName)

        print("This conversion worked both sides.")
        print("Relion values: %s, %s, %s  %s, %s, %s" % (xShift, yShift, zShift, rot, tilt, psi))
        print("Scipion transformation matrix:\n%s\n\n" % np.array_str(alignmentMatrix, precision=2, suppress_small=True))




