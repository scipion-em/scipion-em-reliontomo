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
from os.path import join

import numpy as np
from emtable import Table

from pwem.convert.transformations import translation_from_matrix, euler_from_matrix
from pwem.emlib.image import ImageHandler
from pyworkflow.utils import getExt, removeBaseExt, replaceBaseExt, makePath
from relion.convert.convert_base import WriterBase
from reliontomo import Plugin
from reliontomo.constants import MRC, SHIFTX_ANGST, SHIFTY_ANGST, SHIFTZ_ANGST, TILT, PSI, ROT
from reliontomo.utils import getTransformMatrix
from tomo.objects import Coordinate3D


class WriterTomo(WriterBase):
    def __init__(self,  **kwargs):
        super().__init__(**kwargs)
        self.isPyseg = kwargs.get('isPyseg', False)
        self.starHeaders = kwargs.get('starHeaders', None)
        self.isRelion4 = Plugin.isRe40()


class ReaderTomo:
    def __init__(self, starFile):
        self.starFile = starFile
        self.dataTable = Table()

    def read(self, tableName=None):
        self.dataTable.read(self.starFile, tableName=tableName)


def getTransformInfoFromCoordOrSubtomo(obj, calcInv=True):
    M = obj.getMatrix() if type(obj) is Coordinate3D else obj.getTransform().getMatrix()
    shifts = translation_from_matrix(M)

    if calcInv:
        # Rotation matrix. Remove translation from the Scipion matrix
        R = np.eye(4)
        R[:3, :3] = M[:3, :3]

        Mi = np.linalg.inv(M)

        M = Mi @ R @ R
        shifts = -translation_from_matrix(M)

        M = np.linalg.inv(M)

    angles = -np.rad2deg(euler_from_matrix(M, axes='szyz'))

    return angles, shifts


def checkSubtomogramFormat(subtomo, extraPath):
    """Convert the subtomograms into mrc format if they are in a different one. They will be created in extra folder
    generating a subdirectory for each tomogram found and the corresponding subtomograms contained on it."""
    ih = ImageHandler()
    if getExt(subtomo.getFileName().replace(':' + MRC, '')) != '.' + MRC:
        mrcDir = join(extraPath, removeBaseExt(subtomo.getVolName()))
        makePath(mrcDir)
        mrcFile = join(mrcDir, replaceBaseExt(subtomo.getFileName(), MRC))
        ih.convert(subtomo.getFileName(), mrcFile)


def getTransformMatrixFromRow(row, sRate=1, invert=True):
    shiftx = float(row.get(SHIFTX_ANGST, 0)) / sRate
    shifty = float(row.get(SHIFTY_ANGST, 0)) / sRate
    shiftz = float(row.get(SHIFTZ_ANGST, 0)) / sRate
    tilt = row.get(TILT, 0)
    psi = row.get(PSI, 0)
    rot = row.get(ROT, 0)

    return getTransformMatrix(shiftx, shifty, shiftz, rot, tilt, psi, invert)