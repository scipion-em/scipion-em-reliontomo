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
from enum import Enum
from os import mkdir
from os.path import exists
from pyworkflow.utils import getParentFolder
from reliontomo.constants import SUBTOMO_NAME, FILE_NOT_FOUND, PSUBTOMOS_SQLITE
from reliontomo.convert import createReaderTomo
from reliontomo.objects import RelionSetOfPseudoSubtomograms
from reliontomo.protocols.protocol_base_import_from_star import ProtBaseImportFromStar
from reliontomo.utils import getAbsPath
from tomo.objects import SetOfSubTomograms, SetOfCoordinates3D

SUBTOMOGRAMS_DIR = 'Subtomograms'

class outputObjects(Enum):
    coordinates = SetOfCoordinates3D
    subtomograms = SetOfSubTomograms


class ProtImportSubtomogramsFromStar(ProtBaseImportFromStar):
    """Protocol to import a set of subtomograms from a star file"""

    _label = 'import subtomograms from a star file'
    _possibleOutputs = outputObjects

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    # --------------------------- STEPS functions -----------------------------

    def _initialize(self):
        super()._initialize()
        # Generate the firectoy in which the linked subtomograms pointed from the star file will be stored
        mkdir(self._getExtraPath(SUBTOMOGRAMS_DIR))

    def _importStep(self):
        precedentsSet = self.inTomos.get()
        # Generate the corresponding precedents and 3d coordinates
        super()._importStep()
        # Generate the set of subtomograms
        acq = precedentsSet.getAcquisition()
        sRate = precedentsSet.getSamplingRate()
        if self.readerVersion < 5:
            subtomoSet = SetOfSubTomograms.create(self._getPath(), template='setOfSubTomograms%s.sqlite')
            subtomoSet.setSamplingRate(sRate)
            subtomoSet.setAcquisition(acq)
            self.reader.starFile2SubtomogramsImport(subtomoSet,
                                                    getattr(self, outputObjects.coordinates.name),
                                                    self._getExtraPath(SUBTOMOGRAMS_DIR),
                                                    getParentFolder(self.starFile.get()))
        else:
            subtomoSet = RelionSetOfPseudoSubtomograms.create(self.getPath(), template=PSUBTOMOS_SQLITE)
            subtomoSet.setSamplingRate(sRate)
            subtomoSet.setAcquisition(acq)
            self.reader.starFile2PseudoSubtomograms(subtomoSet)
        self._defineOutputs(**{outputObjects.subtomograms.name: subtomoSet})
        self._defineSourceRelation(self.inTomos, subtomoSet)

    # --------------------------- INFO functions -----------------------------
    def _validate(self):
        errorMsg = super()._validate()
        reader, isReader40 = createReaderTomo(self.starFile.get())
        errorsInPointedFiles = self._checkFilesPointedFromStarFile(getParentFolder(self.starFile.get()),
                                                                   reader.dataTable)
        if errorsInPointedFiles:
            errorMsg.append(errorsInPointedFiles)
        return errorMsg

    # --------------------------- UTILS functions ------------------------------
    def _checkFilesPointedFromStarFile(self, starFilePath, dataTable):
        errorsFound = ''
        fields2check = [SUBTOMO_NAME]
        filesPattern = '\tRow %i - %s\n'

        errorsFound += self._checkFieldsInDataTable(dataTable, fields2check)
        # Check if the files pointed from those fields exist
        if errorsFound:
            filesErrorMsgHead = 'The following files were not found [row, tomoFile, subtomoFile]:\n'
            for counter, row in enumerate(dataTable):
                subtomoFileNotFound = self._fileNotFound(row, SUBTOMO_NAME, starFilePath)
                if subtomoFileNotFound:
                    errorsFound += filesPattern % (counter, subtomoFileNotFound)

            if errorsFound:
                errorsFound = filesErrorMsgHead + errorsFound

        return errorsFound

    @staticmethod
    def _checkFieldsInDataTable(dataTable, fieldList):
        fieldErrors = ''
        if fieldList:
            fieldNotFoundPattern = 'Fields %s were not found in the star file introduced.\n'
            notFoundFields = [field for field in fieldList if not dataTable.hasColumn(field)]
            if notFoundFields:
                pattern = '[%s]' % (' '.join(notFoundFields))
                fieldErrors = (fieldNotFoundPattern % pattern)

        return fieldErrors

    @staticmethod
    def _fileNotFound(row, field, starFilePath):
        statusMsg = ''
        tomoFile = row.get(field, FILE_NOT_FOUND)
        tomoFileAbs = getAbsPath(starFilePath, tomoFile)
        if not exists(tomoFileAbs):
            statusMsg = tomoFile

        return statusMsg
