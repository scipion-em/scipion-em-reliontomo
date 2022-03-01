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
from os import symlink, mkdir
from os.path import exists, join, basename

from pwem.emlib.image import ImageHandler
from pwem.objects import Transform
from pwem.protocols import EMProtocol
from pyworkflow import BETA
from pyworkflow.protocol import FileParam, FloatParam, IntParam, PointerParam
from pyworkflow.utils import Message, removeBaseExt, getParentFolder
from reliontomo.constants import TOMO_NAME_30, PIXEL_SIZE, SUBTOMO_NAME, FILE_NOT_FOUND, TOMO_NAME
from reliontomo.convert import createReaderTomo
from reliontomo.utils import _getAbsPath, manageDims
from tomo.protocols.protocol_base import ProtTomoBase
from tomo.objects import SetOfTomograms, TomoAcquisition, Tomogram, SetOfCoordinates3D


class ProtBaseImportFromStar(EMProtocol, ProtTomoBase):
    """Base protocol for importing data from a star file"""

    _devStatus = BETA

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.linkedStarFileName = None
        self.linkedTomosDirName = 'tomograms'
        self.reader = None
        self.isReader40 = None
        self.sRate = None # Sampling rate of the tomograms where relion picked
        self.starFilePath = None
        self.isSubtomoStarFile = False

    def _defineParams(self, form):
        form.addSection(label=Message.LABEL_INPUT)
        form.addParam('starFile', FileParam,
                      label='Star file')
        form.addParam('inTiltSeries', PointerParam,
                      pointerClass='SetOfTomograms',
                      label='Set of tilt series (opt.)',
                      allowsNull=True,
                      help='Only required if the coordinates are desired to be referred to the corresponding tilt, '
                           'series, like in the case of per-particle per-tilt procedure.')
        form.addParam('samplingRate', FloatParam,
                      label='Sampling rate [Å/pix] (opt.)',
                      allowsNull=True,
                      help='If empty, the protocol will try to read it from the label "_r%s" '
                           'if it is present in the introduced star file.' % PIXEL_SIZE)
        form.addParam('boxSize', IntParam,
                      label='Box Size [pix]',
                      default=20)

    # --------------------------- STEPS functions -----------------------------

    def _insertAllSteps(self):
        self._initialize()
        self._insertFunctionStep(self._importStep)

    def _initialize(self):
        # Get the star file and make a symbolic link in extra directory
        self.starFilePath = getParentFolder(self.starFile.get())
        newStarName = self._getExtraPath(self.linkedStarFileName)
        symlink(self.starFile.get(), newStarName)
        # Read the star file
        self.reader, self.isReader40 = createReaderTomo(starFile=newStarName)
        # Generate the directory in which the linked tomograms pointed from the star file will be stored
        mkdir(self._getExtraPath(self.linkedTomosDirName))
        if self.samplingRate.get():
            self.sRate = self.samplingRate.get()
        else:
            self.sRate = float(self.reader.dataTable[0].get(PIXEL_SIZE))

    def _importStep(self):
        coordSet = SetOfCoordinates3D.create(self._getPath(), template='coordinates%s.sqlite')
        # Generate the precedents (set of tomograms which the coordinates are referred to) if necessary
        if self.isReader40:
            precedentsSet = self.inTiltSeries.get()
            coordSet.setPrecedents(precedentsSet)
        else:
            precedentsSet = SetOfTomograms.create(self._getPath(), template='tomograms%s.sqlite')
            precedentsSet.setSamplingRate(self.sRate)
            precedentsSet.setAcquisition(TomoAcquisition(angleMin=-60, angleMax=60, step=3))  # Generic values
            self._fillPrecedentsSet(precedentsSet)
            self._defineOutputs(outputTomograms=precedentsSet)

        # Read the star file and generate the corresponding set
        # coordSet = self._createSetOfCoordinates3D(precedentsSet)
        coordSet.setSamplingRate(precedentsSet.getSamplingRate())
        coordSet.setBoxSize(self.boxSize.get())
        self.reader.starFile2Coords3D(coordSet, precedentsSet, self.sRate)

        self._defineOutputs(outputCoordinates=coordSet)
        self._defineSourceRelation(precedentsSet, coordSet)

    # --------------------------- INFO functions ------------------------------

    def _validate(self):
        # TODO: check the tsid and tomoname for reader40, CASE PARTIAL MATCHING --> VALIDATE OR WARNING IN SUMMARY
        errors = []
        if not exists(self.starFile.get()):
            errors.append('It was not possible to locate the introduced file. Please check the path.')

        # Check if the files referred in the star file exists
        reader, isReader40 = createReaderTomo(starFile=self.starFile.get())
        errorsInPointedFiles = self._checkFilesPointedFromStarFile(getParentFolder(self.starFile.get()),
                                                                   reader.dataTable,
                                                                   self.isSubtomoStarFile,
                                                                   isReader40)
        if errorsInPointedFiles:
            errors.append(errorsInPointedFiles)

        # In the case of a reader of type 40, the tomoName and the tilt series id must match
        if isReader40:
            if self.inTiltSeries.get():
                tsIds = [ts.getTsId() for ts in self.inTiltSeries.get()]
                if tsIds:
                    coordTomoIds = list(set([row.get(TOMO_NAME) for row in reader.dataTable]))
                    nonMatchingCoorTomoIds = [coordTomoId for coordTomoId in coordTomoIds if coordTomoId not in tsIds]
                    if len(nonMatchingCoorTomoIds) == len(coordTomoIds):
                        errors.append('No matchings were found between the\n'
                                      '\t-Coordinates rlnTomoName [%s]\n'
                                      '\t-Tilt series tsId [%s]' % (', '.join(coordTomoIds), ', '.join(tsIds)))
                else:
                    errors.append('TsId is empty in the introduced tilt series. It is not possible to match the '
                                  'coordintates and the tilt series in that case.')
            else:
                errors.append('A valid set of tilt series is required to match the coordinates in the case of a star '
                              'file of type Relion4.')

        # The sampling rate will be read from the star file if the corresponding field is present in case
        # the user didn't introduce a value for that parameter in the protocol form.
        if not self.samplingRate.get() and not reader.dataTable.hasColumn(PIXEL_SIZE):
            errors.append('Sampling rate was not introduced and the corresponding field [_%s] is not present '
                          'in the introduced star file.' % PIXEL_SIZE)

        return errors

    def _summary(self):
        summary = []
        if self.isFinished():
            if hasattr(self, 'outputTomograms'):
                summary.append('The output set of tomogrmas was generated using the data read from the star file.')
            if not self.samplingRate.get():
                summary.append('The sampling rate considered was the one read from the star file.')

        return summary

    # --------------------------- UTILS functions ------------------------------
    def _fillPrecedentsSet(self, tomoSet):
        """Generate the set of tomograms which the coordinates belong to."""

        tomoFileList = []
        for counter, row in enumerate(self.reader.dataTable):
            tomoFile = row.get(TOMO_NAME_30)
            tomoFileAbs = _getAbsPath(self.starFilePath, tomoFile)
            if exists(tomoFileAbs):
                tomoFileList.append(tomoFile)

        if tomoFileList:
            imgh = ImageHandler()
            # Get the unique elements from the precedent list
            tomoFileList = list(set(tomoFileList))

            for counter, tomoFile in enumerate(tomoFileList):
                tomo = Tomogram()
                origin = Transform()

                # Make a link of the file in extra folder
                linkedFileName = self._getExtraPath(join(self.linkedTomosDirName, basename(tomoFile)))
                symlink(_getAbsPath(self.starFilePath, tomoFile), linkedFileName)

                # Set the tomogram origin
                x, y, z, n = imgh.getDimensions(linkedFileName)
                z = manageDims(linkedFileName, z, n)
                origin.setShifts(x / -2. * self.sRate,
                                 y / -2. * self.sRate,
                                 z / -2. * self.sRate)
                tomo.setOrigin(origin)

                # Set the rest of the tomogram attributes
                tomoId = removeBaseExt(linkedFileName)
                tomo.setTsId(tomoId)
                tomo.setLocation(counter + 1, linkedFileName)
                tomoSet.append(tomo)

        tomoSet.write()
        self._store(tomoSet)

    def _checkFilesPointedFromStarFile(self, starFilePath, dataTable, isSubtomoStarFile=False, isReader40=True):
        errorsFound = ''
        filesPattern = ''
        fields2check = []
        # Check if the corresponding fields exists in the introduced star file
        if isSubtomoStarFile:
            if isReader40:
                filesPattern = '\tRow %i - %s\n'
                fields2check = [SUBTOMO_NAME]
            else:
                filesPattern = '\tRow %i - %s - %s\n'
                fields2check = [TOMO_NAME_30, SUBTOMO_NAME]
        else:
            if not isReader40:
                filesPattern = '\tRow %i - %s\n'
                fields2check = [TOMO_NAME_30]

        errorsFound += self._checkFieldsInDataTable(dataTable, fields2check)
        # Check if the files pointed from those fields exist
        if errorsFound:
            if isSubtomoStarFile:
                filesErrorMsgHead = 'The following files were not found [row, tomoFile, subtomoFile]:\n'
                for counter, row in enumerate(dataTable):
                    subtomoFileNotFound = self._fileNotFound(row, SUBTOMO_NAME, starFilePath)
                    if isReader40:
                        if subtomoFileNotFound:
                            errorsFound += filesPattern % (counter, subtomoFileNotFound)
                    else:
                        tomoFileNotFound = self._fileNotFound(row, TOMO_NAME_30, starFilePath)
                        if tomoFileNotFound or subtomoFileNotFound:
                            errorsFound += filesPattern % (counter, tomoFileNotFound, subtomoFileNotFound)
            else:
                filesErrorMsgHead = 'The following files were not found [row, tomoFile]:\n'
                for counter, row in enumerate(dataTable):
                    if not isReader40:
                        fileNotFound = self._fileNotFound(row, TOMO_NAME_30, starFilePath)
                        if fileNotFound:
                            errorsFound += filesPattern % (counter, fileNotFound)

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
        tomoFileAbs = _getAbsPath(starFilePath, tomoFile)
        if not exists(tomoFileAbs):
            statusMsg = tomoFile

        return statusMsg






