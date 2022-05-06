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
from os import symlink, mkdir
from os.path import exists
from pwem.protocols import EMProtocol
from pyworkflow import BETA
from pyworkflow.protocol import FileParam, FloatParam, IntParam, PointerParam
from pyworkflow.utils import Message, getParentFolder
from reliontomo.constants import PIXEL_SIZE, TOMO_NAME, PARTICLES_TABLE
from reliontomo.convert import createReaderTomo
from tomo.protocols.protocol_base import ProtTomoBase
from tomo.objects import SetOfCoordinates3D


class outputObjects(Enum):
    coordinates = SetOfCoordinates3D


class ProtBaseImportFromStar(EMProtocol, ProtTomoBase):
    """Base protocol for importing data from a star file"""

    _devStatus = BETA

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.linkedStarFileName = None
        self.linkedTomosDirName = 'tomograms'
        self.reader = None
        self.isReader40 = None
        self.coordsSRate = None
        self.starFilePath = None

    def _defineParams(self, form):
        form.addSection(label=Message.LABEL_INPUT)
        form.addParam('starFile', FileParam,
                      label='Star file')
        form.addParam('inTomos', PointerParam,
                      pointerClass='SetOfTomograms',
                      label='Tomograms',
                      help='Tomograms to which the coordinates will be referred.')
        form.addParam('samplingRate', FloatParam,
                      label='Coordinates sampling rate [Å/pix] (opt.)',
                      allowsNull=True,
                      help='If empty, the protocol will try to read it from the label "_r%s" '
                           'if it is present in the introduced star file. If not find, it will be '
                           'considered to be the same as the tomograms. The ratio of both tomograms '
                           'and coordinates sampling rate will be used to scale the coordinates '
                           'properly to the tomograms introduced if they are different' % PIXEL_SIZE)
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
            self.coordsSRate = self.samplingRate.get()
        else:
            self.coordsSRate = float(self.reader.dataTable[0].get(PIXEL_SIZE))

    def _importStep(self):
        precedentsSet = self.inTomos.get()
        coordSet = SetOfCoordinates3D.create(self._getPath(), template='coordinates%s.sqlite')
        coordSet.setPrecedents(precedentsSet)
        coordSet.setSamplingRate(precedentsSet.getSamplingRate())
        coordSet.setBoxSize(self.boxSize.get())
        self.reader.starFile2Coords3D(coordSet, precedentsSet, self.coordsSRate / precedentsSet.getSamplingRate())

        self._defineOutputs(**{outputObjects.coordinates.name: coordSet})
        self._defineSourceRelation(self.inTomos.get(), coordSet)

    # --------------------------- INFO functions ------------------------------

    def _validate(self):
        # TODO: check the tsid and tomoname for reader40, CASE PARTIAL MATCHING --> VALIDATE OR WARNING IN SUMMARY
        errors = []
        if not exists(self.starFile.get()):
            errors.append('It was not possible to locate the introduced file. Please check the path.')

        # Check if the files referred in the star file exists
        reader, isReader40 = createReaderTomo(self.starFile.get())
        if isReader40:
            tsIds = [tomo.getTsId() for tomo in self.inTomos.get()]
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

        return errors

    def _warnings(self):
        warnings = []
        # The sampling rate will be read from the star file if the corresponding field is present in case
        # the user didn't introduce a value for that parameter in the protocol form.
        reader, isReader40 = createReaderTomo(starFile=self.starFile.get())
        if not self.samplingRate.get() and not reader.dataTable.hasColumn(PIXEL_SIZE):
            warnings.append('Sampling rate was not introduced and the corresponding field [_%s] is not present '
                            'in the introduced star file. *Assuming the coordinates to have the same sampling '
                            'rate as the introduced tomograms.*' % PIXEL_SIZE)
        return warnings

    def _summary(self):
        summary = []
        if self.isFinished():
            if hasattr(self, 'outputTomograms'):
                summary.append('The output set of tomogrmas was generated using the data read from the star file.')
            if not self.samplingRate.get():
                summary.append('The sampling rate considered was the one read from the star file.')

        return summary






