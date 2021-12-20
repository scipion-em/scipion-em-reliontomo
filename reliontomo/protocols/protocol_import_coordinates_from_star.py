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
from pyworkflow.protocol import FileParam, FloatParam, IntParam
from pyworkflow.utils import Message, removeBaseExt, getParentFolder
from reliontomo.constants import TOMO_NAME_30, FILE_NOT_FOUND, PIXEL_SIZE
from reliontomo.convert import createReaderTomo30
from reliontomo.utils import _getAbsPath, manageDims
from tomo.protocols.protocol_base import ProtTomoBase
from tomo.objects import SetOfTomograms, TomoAcquisition, Tomogram


class ProtImportCoordinates3DFromStar(EMProtocol, ProtTomoBase):
    """Protocol to import a set of 3D coordinates from a star file"""

    _label = 'import coordinates 3D from a star file'
    _devStatus = BETA

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.linkedStarFileName = 'in3dCoordinates.star'
        self.linkedTomosDirName = 'tomograms'
        self.reader = None
        self.sRate = None
        self.starFilePath = None

    def _defineParams(self, form):
        form.addSection(label=Message.LABEL_INPUT)
        form.addParam('starFile', FileParam,
                      label='3D coordinates star file')
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
        self.reader = createReaderTomo30(starFile=newStarName)
        # Generate the firectoy in which the linked tomograms pointed from the star file will be stored
        mkdir(self._getExtraPath(self.linkedTomosDirName))
        if self.samplingRate.get():
            self.sRate = self.samplingRate.get()
        else:
            self.sRate = float(self.reader.dataTable[0].get(PIXEL_SIZE))

    def _importStep(self):
        # Generate the precedents (set of tomograms which the coordinates are referred to) if necessary
        tomoPrecedentsSet = SetOfTomograms.create(self._getPath(), template='tomograms%s.sqlite')
        tomoPrecedentsSet.setSamplingRate(self.sRate)
        tomoPrecedentsSet.setAcquisition(TomoAcquisition(angleMin=-60, angleMax=60, step=3))  # Generic values
        self._fillPrecedentsSet(tomoPrecedentsSet)
        self._defineOutputs(outputSetOfTomograms=tomoPrecedentsSet)

        # Read the star file and generate the corresponding set
        coordSet = self._createSetOfCoordinates3D(tomoPrecedentsSet)
        coordSet.setSamplingRate(self.sRate)
        coordSet.setBoxSize(self.boxSize.get())
        self.reader.starFile2Coords3D(coordSet, tomoPrecedentsSet)

        self._defineOutputs(outputCoordinates=coordSet)
        self._defineSourceRelation(tomoPrecedentsSet, coordSet)

    # --------------------------- INFO functions ------------------------------

    def _validate(self):
        errors = []
        if not exists(self.starFile.get()):
            errors.append('It was not possible to locate the introduced file. Please check the path.')

        # Check if the files referred in the star file exists
        reader = createReaderTomo30(starFile=self.starFile.get())
        errorsInPointedFiles = self._checkFilesPointedFromStarFile(getParentFolder(self.starFile.get()),
                                                                   reader.dataTable)
        if errorsInPointedFiles:
            errors.append(errorsInPointedFiles)

        # The sampling rate will be read from the star file if the corresponding field is present in case
        # the user didn't introduce a value for that parameter in the protocol form.
        if not self.samplingRate.get() and not reader.dataTable.hasColumn(PIXEL_SIZE):
            errors.append('Sampling rate was not introduced and the corresponding field [_%s] is not present '
                          'in the introduced star file.' % PIXEL_SIZE)

        return errors

    def _summary(self):
        summary = []
        if self.isFinished():
            summary.append('The output set of tomogrmas was generated using the data read from the star file.')
            if not self.samplingRate.get():
                summary.append('The sampling rate considered was the one read from the star file.')

        return summary

    # --------------------------- UTILS functions ------------------------------

    @staticmethod
    def _checkFilesPointedFromStarFile(starFilePath, dataTable):
        tomogramsNotFoundMsg = ''
        fieldNotFoundPattern = 'Field %s was not found in the star file introduced.\n'
        if not dataTable.hasColumn(TOMO_NAME_30):
            tomogramsNotFoundMsg += fieldNotFoundPattern % TOMO_NAME_30
        if dataTable.hasColumn(TOMO_NAME_30):
            OneFilePattern = '\tRow %i - %s\n'
            for counter, row in enumerate(dataTable):
                tomoFile = row.get(TOMO_NAME_30, FILE_NOT_FOUND)
                tomoFileAbs = _getAbsPath(starFilePath, tomoFile)
                if not exists(tomoFileAbs):
                    tomogramsNotFoundMsg += OneFilePattern % (counter, tomoFile)

        if tomogramsNotFoundMsg:
            tomogramsNotFoundMsg = 'The following files were not found [row, tomogram]:\n\n' + \
                                   tomogramsNotFoundMsg
        return tomogramsNotFoundMsg

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






