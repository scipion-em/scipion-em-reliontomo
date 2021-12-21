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


def getProgram(program, nMpi):
    """ Get the program name depending on the MPI use or not."""
    if nMpi > 1:
        program += '_mpi'
    return program


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


# def _checkFilesPointedFromStarFile(starFilePath, dataTable, isSubtomoStarFile=False):
#     from reliontomo.constants import TOMO_NAME_30, SUBTOMO_NAME
#     errorsFound = ''
#     # Check if the corresponding fields exists in the introduced star file
#     fields2check = [TOMO_NAME_30]
#     filesPattern = '\tRow %i - %s\n'
#     if isSubtomoStarFile:
#         fields2check = [TOMO_NAME_30, SUBTOMO_NAME]
#         filesPattern = '\tRow %i - %s - %s\n'
#
#     errorsFound += _checkFieldsInDataTable(dataTable, fields2check)
#     # Check if the files pointed from those fields exist
#     if not errorsFound:
#         if isSubtomoStarFile:
#             filesErrorMsgHead = 'The following files were not found [row, tomoFile, subtomoFile]:\n'
#             for counter, row in enumerate(dataTable):
#                 tomoFileNotFound = _fileNotFound(row, TOMO_NAME_30, starFilePath)
#                 subtomoFileNotFound = _fileNotFound(row, SUBTOMO_NAME, starFilePath)
#                 if tomoFileNotFound or subtomoFileNotFound:
#                     errorsFound += filesPattern % (counter, tomoFileNotFound. subtomoFileNotFound)
#         else:
#             filesErrorMsgHead = 'The following files were not found [row, tomoFile]:\n'
#             for counter, row in enumerate(dataTable):
#                 fileNotFound = _fileNotFound(row, TOMO_NAME_30, starFilePath)
#                 if fileNotFound:
#                     errorsFound += filesPattern % (counter, fileNotFound)
#
#         if errorsFound:
#             errorsFound = filesErrorMsgHead + errorsFound
#
#     return errorsFound
#
#
# def _checkFieldsInDataTable(dataTable, fieldList):
#     fieldErrors = ''
#     fieldNotFoundPattern = 'Fields %s were not found in the star file introduced.\n'
#     notFoundFields = [field for field in fieldList if not dataTable.hasColumn(field)]
#     if notFoundFields:
#         pattern = '[%s]' % (' '.join(notFoundFields))
#         fieldErrors = (fieldNotFoundPattern % pattern)
#
#     return fieldErrors
#
#
# def _fileNotFound(row, field, starFilePath):
#     from reliontomo.constants import FILE_NOT_FOUND
#     statusMsg = ''
#     tomoFile = row.get(field, FILE_NOT_FOUND)
#     tomoFileAbs = _getAbsPath(starFilePath, tomoFile)
#     if not exists(tomoFileAbs):
#         statusMsg = tomoFile
#
#     return statusMsg

