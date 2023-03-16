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
from os.path import exists, join

from pwem.objects import VolumeMask
from pwem.protocols import EMProtocol
from pyworkflow.protocol import PointerParam
from pyworkflow.utils import Message, createLink
from reliontomo.constants import IN_PARTICLES_STAR, POSTPROCESS_DIR, OPTIMISATION_SET_STAR, PSUBTOMOS_SQLITE
from reliontomo.convert import writeSetOfPseudoSubtomograms, readSetOfPseudoSubtomograms
from reliontomo.objects import createSetOfRelionPSubtomograms, RelionSetOfPseudoSubtomograms


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

    def genInStarFile(self, withPriors=False):
        """It will check if the set size and the stored particles star file are of the same size or not. In
        the first case, a link will be made to the previous particles star file to avoid generating it and in the
        second case, a new file will be generated containing only the ones present in the input set."""
        inReParticlesSet = self.inReParticles.get()
        outStarFileName = self.getOutStarFileName()
        # if inReParticlesSet.getSize() == inReParticlesSet.getNReParticles():
        #     createLink(inReParticlesSet.getParticles(), outStarFileName)
        # else:
        writeSetOfPseudoSubtomograms(inReParticlesSet, outStarFileName, withPriors=withPriors)

    def _genPostProcessOutputMrcFile(self, fileName):
        """File generated using the sharpening protocol (called post-process protocol) and also using the
        rec particle from TS protocol in case the optional input 'solvent mask' is introduced."""
        postProccesMrc = VolumeMask()
        postProccesMrc.setFileName(self._getExtraPath(POSTPROCESS_DIR, fileName))
        postProccesMrc.setSamplingRate(self.inReParticles.get().getCurrentSamplingRate())

        return postProccesMrc

    def genRelionParticles(self, binningFactor=None, boxSize=24):
        """Generate a RelionSetOfPseudoSubtomograms object containing the files involved for the next protocol,
        considering that some protocols don't generate the optimisation_set.star file. In that case, the input Object
        which represents it will be copied and, after that, this method will be used to update the corresponding
        attribute."""
        extraPath = self._getExtraPath()
        inParticlesSet = self.inReParticles.get()
        coordsPointer = inParticlesSet._coordsPointer
        optimSetStar = join(extraPath, OPTIMISATION_SET_STAR)
        protocolPath = join(extraPath, '..')
        if exists(optimSetStar):
            psubtomoSet = createSetOfRelionPSubtomograms(protocolPath,
                                                         optimSetStar,
                                                         coordsPointer,
                                                         template=PSUBTOMOS_SQLITE,
                                                         tsSamplingRate=inParticlesSet.getTsSamplingRate(),
                                                         relionBinning=binningFactor if binningFactor else inParticlesSet.getRelionBinning(),
                                                         boxSize=boxSize if boxSize else inParticlesSet.getBoxSize())
        else:
            psubtomoSet = RelionSetOfPseudoSubtomograms.create(protocolPath, template=PSUBTOMOS_SQLITE)
            psubtomoSet.copyInfo(inParticlesSet)
            psubtomoSet.updateGenFiles(extraPath)

        # Fill the set with the generated particles
        readSetOfPseudoSubtomograms(psubtomoSet)
        return psubtomoSet
