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
import json
from os.path import isabs, abspath, join


def getProgram(program, nMpi):
    """ Get the program name depending on the MPI use or not."""
    if nMpi > 1:
        program += '_mpi'
    return program


def genSymmetryTable():
    jsonData = '[{"group": "Asymmetric", "notation": "C1", "origin": "User-defined", "orientation": "User-defined"},' \
               '{"group": "Cyclic", "notation": "C<n>", "origin": "On symm axis, Z user-defined", "orientation": "Symm axis on Z"},' \
               '{"group": "Dihedral", "notation": "D<n>", "origin": "Intersection of symm axes", "orientation": "principle symm axis on Z, 2-fold on X"},' \
               '{"group": "Tetrahedral", "notation": "T", "origin": "Intersection of symm axes", "orientation": "3-fold axis on Z (deviating from Heymann et al!)"},' \
               '{"group": "Octahedral", "notation": "O", "origin": "Intersection of symm axes", "orientation": "4-fold axes on X, Y, Z"},' \
               '{"group": "Icosahedral", "notation": "I<n>", "origin": "Intersection of symm axes", "orientation": "**"}]'

    return json.loads(jsonData)


def getFileFromDataPrepProt(prot, fileName):
    return prot.inputPrepareDataProt.get()._getExtraPath(fileName)


def manageDims(fileName, z, n):
    if fileName.endswith('.mrc') or fileName.endswith('.map'):
        if z == 1 and n != 1:
            zDim = n
        else:
            zDim = z
    else:
        zDim = z

    return zDim


def _getAbsPath(starFilePath, tomoFile):
    """If the paths of the files pointed from a star file are relative, they'll be referred to the path of the
    star file. This method is used to consider that case."""
    if isabs(tomoFile):
        return tomoFile
    else:
        return join(starFilePath, tomoFile)

