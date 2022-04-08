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
from os.path import isabs, join, exists
import numpy as np
from pwem.convert import transformations
from reliontomo.constants import OPTIMISATION_SET_STAR, OUT_PARTICLES_STAR, PSUBTOMOS_SQLITE
from reliontomo.objects import relionTomoMetadata, SetOfPseudoSubtomograms


def getProgram(program, nMpi):
    """ Get the program name depending on the MPI use or not."""
    if nMpi > 1:
        program += '_mpi'
    return program


def manageDims(fileName, z, n):
    if fileName.endswith('.mrc') or fileName.endswith('.map'):
        if z == 1 and n != 1:
            zDim = n
        else:
            zDim = z
    else:
        zDim = z

    return zDim


def getAbsPath(starFilePath, tomoFile):
    """If the paths of the files pointed from a star file are relative, they'll be referred to the path of the
    star file. This method is used to consider that case."""
    if isabs(tomoFile):
        return tomoFile
    else:
        return join(starFilePath, tomoFile)


def _gen2LevelBaseName(fullFileName):
    """This function generates a new basename composed of the original basename preceded by the parent folder
    followed by '_'. This will reduce the probabilities of non-uniqueness of the name when creating links."""
    parts = fullFileName.split('/')
    if len(parts) == 1:
        return fullFileName
    else:
        return parts[-2] + '_' + parts[-1]


def isPseudoSubtomogram(subtomo):
    return hasattr(subtomo, '_ctfImage')


def getTransformMatrix(shiftx, shifty, shiftz, rot, tilt, psi, invert):
    shifts = (float(shiftx), float(shifty), float(shiftz))
    angles = (float(rot), float(tilt), float(psi))
    radAngles = -np.deg2rad(angles)
    M = transformations.euler_matrix(radAngles[0], radAngles[1], radAngles[2], 'szyz')
    if invert:
        M[0, 3] = -shifts[0]
        M[1, 3] = -shifts[1]
        M[2, 3] = -shifts[2]
        M = np.linalg.inv(M)
    else:
        M[0, 3] = shifts[0]
        M[1, 3] = shifts[1]
        M[2, 3] = shifts[2]

    return M


def genRelionParticles(extraPath, inOptSet, binningFactor=None, nParticles=None):
    """Generate a relionParticles object containing the files involved for the next protocol, considering that some
    protocols don't generate the optimisation_set.star file. In that case, the input Object which represents it will
    be copied and, after that, this method will be used to update the corresponding attribute."""
    optimSetStar = join(extraPath, OPTIMISATION_SET_STAR)
    if exists(optimSetStar):
        relionParticles = relionTomoMetadata(optimSetStar=optimSetStar,
                                             tsSamplingRate=inOptSet.getTsSamplingRate(),
                                             relionBinning=binningFactor if binningFactor else inOptSet.getRelionBinning(),
                                             nParticles=nParticles if nParticles else inOptSet.getNumParticles())
    else:
        relionParticles = relionTomoMetadata()
        relionParticles.copyInfo(inOptSet)
        relionParticles.updateGenFiles(extraPath)

    return relionParticles


def genOutputPseudoSubtomograms(prot, currentSamplingRate):
    """Centralized code to generate the output set of pseudosubtomograms for protocols make pseudosubtomograms,
     auto-refine, CTF refine and frame align"""
    from reliontomo.convert import createReaderTomo
    reader, _ = createReaderTomo(prot._getExtraPath(OUT_PARTICLES_STAR))
    outputSet = prot._createSet(SetOfPseudoSubtomograms, PSUBTOMOS_SQLITE, '')
    outputSet.setSamplingRate(currentSamplingRate)
    reader.starFile2PseudoSubtomograms(outputSet)
    return outputSet
