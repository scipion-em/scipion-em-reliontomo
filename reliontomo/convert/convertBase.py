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
import logging
from os.path import join
import numpy as np
from pwem.convert import transformations
from pwem.convert.transformations import translation_from_matrix, euler_from_matrix
from pwem.emlib.image import ImageHandler
from pyworkflow.utils import getExt, removeBaseExt, replaceBaseExt, makePath, cyanStr
from relion.convert.convert_base import WriterBase
from reliontomo.constants import MRC, SHIFTX_ANGST, SHIFTY_ANGST, SHIFTZ_ANGST, TILT, PSI, ROT
from tomo.constants import TR_RELION
from tomo.objects import Coordinate3D

logger = logging.getLogger(__name__)


class WriterTomo(WriterBase):
    def __init__(self,  **kwargs):
        super().__init__(**kwargs)
        self.starHeaders = kwargs.get('starHeaders', None)


class ReaderTomo:
    def __init__(self, starFile, dataTable):
        self.starFile = starFile
        self.dataTable = dataTable


def getTransformInfoFromCoordOrSubtomo(obj, samplingRate):

    M = obj.getMatrix(convention=TR_RELION) if isinstance(obj, Coordinate3D) else obj.getTransform(convention=TR_RELION).getMatrix()
    shifts = translation_from_matrix(M)

    # These 2 lines below were done when inverting, which is now what we always do.
    shifts = -shifts
    M = np.linalg.inv(M)

    angles = -np.rad2deg(euler_from_matrix(M, axes='szyz'))
    shifts *= samplingRate

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


# def getTransformMatrixFromRow(row, sRate=1, isRe5Star=False):
#     if isRe5Star:
#         from reliontomo.convert.convert50_tomo import RLN_ORIGINZANGST, RLN_ORIGINYANGST, RLN_ORIGINXANGST, \
#             RLN_TOMOSUBTOMOGRAMROT, RLN_TOMOSUBTOMOGRAMTILT, RLN_TOMOSUBTOMOGRAMPSI
#         shiftx = float(row.get(RLN_ORIGINXANGST, 0))
#         shifty = float(row.get(RLN_ORIGINYANGST, 0))
#         shiftz = float(row.get(RLN_ORIGINZANGST, 0))
#         rot = row.get(RLN_TOMOSUBTOMOGRAMROT, 0)
#         tilt = row.get(RLN_TOMOSUBTOMOGRAMTILT, 0)
#         psi = row.get(RLN_TOMOSUBTOMOGRAMPSI, 0)
#     else:
#         shiftx = float(row.get(SHIFTX_ANGST, 0))
#         shifty = float(row.get(SHIFTY_ANGST, 0))
#         shiftz = float(row.get(SHIFTZ_ANGST, 0))
#         rot = row.get(ROT, 0)
#         tilt = row.get(TILT, 0)
#         psi = row.get(PSI, 0)
#
#     return genTransformMatrix(shiftx, shifty, shiftz, rot, tilt, psi, sRate)


def getTransformMatrixFromRow(row, sRate=1, isRe5Star=False):
    if isRe5Star:
        logger.info(cyanStr('Is Relion 5'))
        from reliontomo.convert.convert50_tomo import RLN_ORIGINZANGST, RLN_ORIGINYANGST, RLN_ORIGINXANGST, \
            RLN_TOMOSUBTOMOGRAMROT, RLN_TOMOSUBTOMOGRAMTILT, RLN_TOMOSUBTOMOGRAMPSI
        shiftx = float(row.get(RLN_ORIGINXANGST, 0))
        shifty = float(row.get(RLN_ORIGINYANGST, 0))
        shiftz = float(row.get(RLN_ORIGINZANGST, 0))
        rot = row.get(RLN_TOMOSUBTOMOGRAMROT, 0)
        tilt = row.get(RLN_TOMOSUBTOMOGRAMTILT, 0)
        psi = row.get(RLN_TOMOSUBTOMOGRAMPSI, 0)
    else:
        logger.info(cyanStr('It is NOT Relion 5'))
        shiftx = float(row.get(SHIFTX_ANGST, 0))
        shifty = float(row.get(SHIFTY_ANGST, 0))
        shiftz = float(row.get(SHIFTZ_ANGST, 0))
        rot = row.get(ROT, 0)
        tilt = row.get(TILT, 0)
        psi = row.get(PSI, 0)

    logger.info(cyanStr(f'rot = {rot}, tilt = {tilt}, psi = {psi}'))
    return genTransformMatrix(shiftx, shifty, shiftz, rot, tilt, psi, sRate)


def genTransformMatrix(shiftx, shifty, shiftz, rot, tilt, psi, sRate):
    shifts = (float(shiftx)/sRate, float(shifty)/sRate, float(shiftz)/sRate)
    angles = (float(rot), float(tilt), float(psi))
    radAngles = -np.deg2rad(angles)
    M = transformations.euler_matrix(radAngles[0], radAngles[1], radAngles[2], 'szyz')

    # These 3 lines are the ones for "invert" flag.
    M[0, 3] = -shifts[0]
    M[1, 3] = -shifts[1]
    M[2, 3] = -shifts[2]
    M = np.linalg.inv(M)

    # These line were fot the non invert mode. Not used
    #     M[0, 3] = shifts[0]
    #     M[1, 3] = shifts[1]
    #     M[2, 3] = shifts[2]

    return M
