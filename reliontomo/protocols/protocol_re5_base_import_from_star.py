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
from os.path import exists
from pwem.protocols import EMProtocol
from pyworkflow.protocol import FileParam, FloatParam, IntParam, PointerParam
from pyworkflow.utils import Message, getParentFolder, createLink, makePath
from reliontomo.constants import TOMO_NAME
from reliontomo.convert import createReaderTomo
from tomo.protocols.protocol_base import ProtTomoBase
from tomo.objects import SetOfCoordinates3D


class outputObjects(Enum):
    coordinates = SetOfCoordinates3D


class ProtBaseRe5ImportFromStar(EMProtocol, ProtTomoBase):
    """Base protocol for importing data from a star file"""

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.linkedStarFileName = None
        self.linkedTomosDirName = 'tomograms'
        self.reader = None
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
                      help='If empty, it will be considered to be the same as the tomograms. The ratio of both '
                           'tomograms and coordinates sampling rate will be used to scale the coordinates '
                           'properly to the tomograms introduced if they are different.')
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
        createLink(self.starFile.get(), newStarName)
        # Read the star file
        self.reader, _ = createReaderTomo(starFile=newStarName)
        # Generate the directory in which the linked tomograms pointed from the star file will be stored
        makePath(self._getExtraPath(self.linkedTomosDirName))
        if self.samplingRate.get():
            self.coordsSRate = self.samplingRate.get()
        else:
            self.coordsSRate = self.inTomos.get().getSamplingRate()

    def _importStep(self):
        inTomosPointer = self.inTomos
        precedentsSet = inTomosPointer.get()
        coordSet = SetOfCoordinates3D.create(self._getPath(), template='coordinates%s.sqlite')
        coordSet.setPrecedents(inTomosPointer)  # As a pointer is better for streaming
        coordSet.setSamplingRate(precedentsSet.getSamplingRate())
        coordSet.setBoxSize(self.boxSize.get())
        self.reader.starFile2Coords3D(coordSet, precedentsSet, self.coordsSRate / precedentsSet.getSamplingRate())

        self._defineOutputs(**{outputObjects.coordinates.name: coordSet})
        self._defineSourceRelation(inTomosPointer, coordSet)

    # --------------------------- INFO functions ------------------------------

    def _validate(self):
        errors = []
        if not exists(self.starFile.get()):
            errors.append('It was not possible to locate the introduced file. Please check the path.')

        # Check if the files referred in the star file exists
        reader, isReaderGE40 = createReaderTomo(self.starFile.get())
        if isReaderGE40:
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

    def _summary(self):
        summary = []
        if self.isFinished():
            if not self.samplingRate.get():
                summary.append('The sampling rate considered was the one read from the star file.')
            summary.append("Coordinates imported from %s." % self.starFile.get())

        return summary
