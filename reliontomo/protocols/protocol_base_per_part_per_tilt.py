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
from pyworkflow.protocol import PointerParam, IntParam, GE, LE
from reliontomo.objects import RelionSetOfPseudoSubtomograms
from reliontomo.protocols.protocol_base_relion import ProtRelionTomoBase


class outputObjects(Enum):
    relionParticles = RelionSetOfPseudoSubtomograms


class ProtRelionPerParticlePerTiltBase(ProtRelionTomoBase):
    """Base protocol used for the getting the frame alignment and ctf-refinment"""

    _possibleOutputs = outputObjects

    # -------------------------- DEFINE param functions -----------------------

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.inParticlesStar = None

    def _defineParams(self, form):
        super()._defineCommonInputParams(form)
        form.addParam('recVolume', PointerParam,
                      pointerClass='AverageSubTomogram',
                      allowsNull=False,
                      label='Volume to get the halves')
        form.addParam('inRefMask', PointerParam,
                      pointerClass='VolumeMask',
                      label="Input reference mask")

    @staticmethod
    def _insertBoxSizeForEstimationParam(form):
        form.addParam('boxSize', IntParam,
                      label='Box size for estimation (pix)',
                      default=128,
                      allowsNull=False,
                      validators=[GE(32), LE(512)],
                      help="Box size to be used for the estimation. Note that this can be larger than the box size "
                           "of the reference map. A sufficiently large box size allows more of the high-frequency "
                           "signal to be captured that has been delocalized by the CTF.")

    # -------------------------- INSERT steps functions -----------------------
    def _insertAllSteps(self):
        pass

    def createOutputStep(self):
        inPSubtomos = self.inReParticles.get()
        # Output RelionParticles
        pSubtomos = self.genRelionParticles()

        self._defineOutputs(**{outputObjects.relionParticles.name: pSubtomos})
        self._defineSourceRelation(inPSubtomos, pSubtomos)

    # -------------------------- UTILS functions -----------------------------
    def _genIOCommand(self):
        inPSubtomos = self.inReParticles.get()
        trajectories = inPSubtomos.getTrajectories()
        postProcess = inPSubtomos.getReferenceFsc()
        half1, half2 = self.recVolume.get().getHalfMaps().split(',')
        cmd = '--p %s ' % inPSubtomos.getParticles()
        cmd += '--t %s ' % inPSubtomos.getTomogramsStar()
        cmd += '--o %s ' % self._getExtraPath()
        if trajectories:
            cmd += '--mot %s ' % trajectories
        cmd += '--ref1 %s ' % half1
        cmd += '--ref2 %s ' % half2
        cmd += '--mask %s ' % self.inRefMask.get().getFileName()
        if postProcess:
            cmd += '--fsc %s ' % postProcess
        cmd += '--b %i ' % self.boxSize.get()
        cmd += self._genExtraParamsCmd()
        return cmd







