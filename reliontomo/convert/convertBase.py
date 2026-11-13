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


# ------------------------------------------------------------------------------------------------------------------
# Relion 5 particles: a particle is described by two sets of angles and one shift vector (see Relion's
# src/jaz/tomography/particle_set.cpp):
#   - A_sub  <-- rlnTomoSubtomogramRot|Tilt|Psi: prior orientation (e.g. oriented picking). Not refined.
#   - A_part <-- rlnAngleRot|Tilt|Psi: refined orientation, relative to A_sub.
#   - o      <-- rlnOriginX|Y|ZAngst: refined shifts, expressed in the frame defined by A_sub.
# so that:
#   particle orientation = A_sub @ A_part
#   particle position    = coordinate - A_sub @ o / pixelSize
# In Scipion (particle(r) = reference(M @ r)) they are stored as:
#   - Coordinate matrix (subtomo.getCoordinate3D().getMatrix()), the prior: [ A_sub^-1 | 0 ]
#   - Transform (subtomo.getTransform()), absolute, the prior and the refinement: [ (A_sub @ A_part)^-1 | A_part^-1 @ o ]
# ------------------------------------------------------------------------------------------------------------------
def relionAngles2Matrix(rot, tilt, psi):
    """Returns the 3x3 rotation matrix A that corresponds to the given Relion euler angles (in degrees)."""
    radAngles = -np.deg2rad((float(rot), float(tilt), float(psi)))
    return transformations.euler_matrix(radAngles[0], radAngles[1], radAngles[2], 'szyz')[:3, :3]


def matrix2RelionAngles(A, positiveTilt=False):
    """Returns the Relion euler angles (rot, tilt, psi, in degrees) that correspond to the 3x3 rotation matrix A.
    :param positiveTilt: if True, the equivalent solution with the tilt angle in the range [0, 180] is returned."""
    M = np.eye(4)
    M[:3, :3] = A
    rot, tilt, psi = -np.rad2deg(euler_from_matrix(M, axes='szyz'))
    if positiveTilt and tilt < 0:
        # (rot, tilt, psi) and (rot + 180, -tilt, psi - 180) represent the same rotation
        rot, tilt, psi = rot + 180, -tilt, psi - 180
        rot = (rot + 180) % 360 - 180
        psi = (psi + 180) % 360 - 180
    return np.array([rot, tilt, psi])


def genRe5ParticleMatrices(subtomoAngles, particleAngles, shiftsAngst, sRate):
    """Generates the Scipion matrices of a Relion 5 particle.
    :param subtomoAngles: rlnTomoSubtomogramRot|Tilt|Psi, in degrees.
    :param particleAngles: rlnAngleRot|Tilt|Psi, in degrees.
    :param shiftsAngst: rlnOriginX|Y|ZAngst, in angstroms.
    :param sRate: sampling rate of the particles, in angstroms/pixel.
    :return: the matrix of the coordinate (prior orientation) and the matrix of the particle (absolute: prior
    orientation and refinement, with the shifts in pixels)."""
    aSub = relionAngles2Matrix(*subtomoAngles)
    aPart = relionAngles2Matrix(*particleAngles)
    shifts = np.array(shiftsAngst, dtype=float) / sRate
    coordMatrix = np.eye(4)
    coordMatrix[:3, :3] = aSub.T
    particleMatrix = np.eye(4)
    particleMatrix[:3, :3] = (aSub @ aPart).T
    particleMatrix[:3, 3] = aPart.T @ shifts
    return coordMatrix, particleMatrix


def getRe5ParticleAnglesAndShifts(coordMatrix, particleMatrix, sRate, refAngles=None, refShiftsAngst=None,
                                  angTol=1e-3, shiftTol=1e-3):
    """Inverse operation of method genRe5ParticleMatrices.
    :param coordMatrix: Scipion matrix of the coordinate (prior orientation).
    :param particleMatrix: Scipion matrix of the particle (absolute).
    :param sRate: sampling rate of the particles, in angstroms/pixel.
    :param refAngles: rlnAngleRot|Tilt|Psi values stored with the particle. They are returned unchanged if they
    represent the rotation contained in the matrices (the euler angles that represent a rotation are not unique),
    so the values written by Relion are preserved unless the matrices were modified.
    :param refShiftsAngst: the same, for rlnOriginX|Y|ZAngst.
    :param angTol: angular tolerance, in degrees, used to compare with refAngles.
    :param shiftTol: tolerance, in angstroms, used to compare with refShiftsAngst.
    :return: subtomoAngles (deg), particleAngles (deg) and shifts (angstroms)."""
    aSub = coordMatrix[:3, :3].T
    aPart = aSub.T @ particleMatrix[:3, :3].T
    shifts = (aPart @ particleMatrix[:3, 3]) * sRate
    subtomoAngles = matrix2RelionAngles(aSub)
    particleAngles = None
    if refAngles is not None and None not in refAngles:
        diffRot = relionAngles2Matrix(*refAngles).T @ aPart
        diffAngle = np.rad2deg(np.arccos(np.clip((np.trace(diffRot) - 1) / 2, -1, 1)))
        if diffAngle <= angTol:
            particleAngles = np.array(refAngles, dtype=float)
    if particleAngles is None:
        particleAngles = matrix2RelionAngles(aPart, positiveTilt=True)
    if refShiftsAngst is not None and None not in refShiftsAngst:
        if np.abs(shifts - np.array(refShiftsAngst, dtype=float)).max() <= shiftTol:
            shifts = np.array(refShiftsAngst, dtype=float)
    return subtomoAngles, particleAngles, shifts


def getTransformMatrixFromRow(row, sRate=1, isRe5Star=False):
    if isRe5Star:
        from reliontomo.convert.convert50_tomo import RLN_ORIGINZANGST, RLN_ORIGINYANGST, RLN_ORIGINXANGST, \
            RLN_TOMOSUBTOMOGRAMROT, RLN_TOMOSUBTOMOGRAMTILT, RLN_TOMOSUBTOMOGRAMPSI
        shiftx = float(row.get(RLN_ORIGINXANGST, 0))
        shifty = float(row.get(RLN_ORIGINYANGST, 0))
        shiftz = float(row.get(RLN_ORIGINZANGST, 0))
        rot = row.get(RLN_TOMOSUBTOMOGRAMROT, 0)
        tilt = row.get(RLN_TOMOSUBTOMOGRAMTILT, 0)
        psi = row.get(RLN_TOMOSUBTOMOGRAMPSI, 0)
    else:
        shiftx = float(row.get(SHIFTX_ANGST, 0))
        shifty = float(row.get(SHIFTY_ANGST, 0))
        shiftz = float(row.get(SHIFTZ_ANGST, 0))
        rot = row.get(ROT, 0)
        tilt = row.get(TILT, 0)
        psi = row.get(PSI, 0)

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
