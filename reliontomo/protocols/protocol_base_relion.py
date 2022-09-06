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
from pwem.protocols import EMProtocol
from pyworkflow.protocol import PointerParam
from pyworkflow.utils import Message, createLink
from reliontomo.constants import IN_PARTICLES_STAR
from reliontomo.convert import writeSetOfPseudoSubtomograms


class ProtRelionTomoBase(EMProtocol):

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    def _defineCommonInputParams(self, form):
        form.addSection(label=Message.LABEL_INPUT)
        form.addParam('inReParticles', PointerParam,
                      pointerClass='RelionSetOfPseudoSubtomograms',
                      label='Relion particles')

    def getOutStarFileName(self):
        return self._getExtraPath(IN_PARTICLES_STAR)

    def genInStarFile(self):
        """It will check if the set size and the stored particles star file are of the same size or not. In
        the first case, a link will be made to the previous particles star file to avoid generating it and in the
        second case, a new file will be generated containing only the ones present in the input set."""
        inReParticlesSet = self.inReParticles.get()
        outStarFileName = self.getOutStarFileName()
        if inReParticlesSet.getSize() == inReParticlesSet.getNReParticles():
            createLink(inReParticlesSet.getParticles(), outStarFileName)
        else:
            writeSetOfPseudoSubtomograms(inReParticlesSet, outStarFileName)
