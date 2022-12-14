# **************************************************************************
# *
# * Authors: Yunior C. Fonseca Reyna    (cfonseca@cnb.csic.es)
# *
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
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************
from enum import Enum
from pyworkflow import BETA
from pyworkflow.protocol.params import (PointerParam)
from pwem.protocols import EMProtocol
from reliontomo.objects import RelionSetOfPseudoSubtomograms

from tomo.protocols import ProtTomoBase
from tomo.objects import SetOfSubTomograms


class outputObjects(Enum):
    relionParticles = RelionSetOfPseudoSubtomograms

class MatchingCoordinates(EMProtocol, ProtTomoBase):
    """
    """
    _label = 'Matching coordinates'
    _devStatus = BETA

    # --------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputPseudoSubtomograms', PointerParam,
                      pointerClass='RelionSetOfPseudoSubtomograms',
                      label="Input pseudoubtomograms", important=True,
                      help="Select the input pseudosubtomograms from the project.")
        form.addParam('inputSubTomos3DCoord', PointerParam,
                      pointerClass='SetOfSubTomograms, SetOfCoordinates3D',
                      label='Subtomograms or 3D coordinates', important=True,
                      help='Select the subtomograms or 3D coordinates that we '
                           'want to match with the pseudosubtomogram')

    # --------------------------- INSERT steps functions ----------------------
    def _insertAllSteps(self):
        self._insertFunctionStep(self.createOutputStep)

    # --------------------------- STEPS functions -----------------------------

    def createOutputStep(self):
        inputPseudoSubtomograms = self.inputPseudoSubtomograms.get()
        pseudoSubTSr = inputPseudoSubtomograms.getSamplingRate()
        inputSubTomosOr3DCoord = self.inputSubTomos3DCoord.get()
        subtomoSr = inputSubTomosOr3DCoord.getSamplingRate()


        if isinstance(inputSubTomosOr3DCoord, SetOfSubTomograms):
            inputSubTomosOr3DCoord = inputSubTomosOr3DCoord.getCoordinates3D()

        protocolPath =self.getPath()
        psubtomoSet = RelionSetOfPseudoSubtomograms.create(protocolPath,
                                                           template="pseudosubtomograms.sqlite")
        psubtomoSet.copyInfo(inputPseudoSubtomograms)

        subtomo3DCoordList = [subtomo3DCoord.composeCoordId(subtomoSr) for subtomo3DCoord in inputSubTomosOr3DCoord.iterCoordinates()]
        seudoSubTomoDict = dict()
        for seudoSubTomo in inputPseudoSubtomograms:
            coord = seudoSubTomo.getCoordinate3D()
            seudoSubTomoDict[coord.composeCoordId(pseudoSubTSr)] = seudoSubTomo.clone()

        for subtomo3DCoord in subtomo3DCoordList:
            if subtomo3DCoord in seudoSubTomoDict:
                psubtomoSet.append(seudoSubTomoDict[subtomo3DCoord])

        self._defineOutputs(**{outputObjects.relionParticles.name: psubtomoSet})





