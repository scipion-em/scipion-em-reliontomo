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
from pyworkflow.protocol import IntParam, FloatParam, GE
from reliontomo.constants import IN_PARTICLES_STAR, IN_TOMOS_STAR
from reliontomo.objects import RelionSetOfPseudoSubtomograms
from reliontomo.protocols.protocol_base_relion import ProtRelionTomoBase


class ProtRelion5ExtractSubtomoAndRecParticleBase(ProtRelionTomoBase):
    """Reconstruct particle and make pseudo-subtomograms base class"""

    _label = None

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        self._defineCommonInputParams(form)
        form.addParallelSection(threads=1, mpi=3)

    @staticmethod
    def _defineCommonRecParams(form):
        form.addParam('binningFactor', FloatParam,
                      label='Binning factor (downsampling)',
                      default=1,
                      validator=GE(0),
                      allowsNull=False,
                      help='The tilt series images will be binned by this (real-valued) factor and then '
                           'reconstructed in the specified box size above. Note that thereby the '
                           'reconstructed region becomes larger when specifying binning factors larger than one. '
                           'This does not alter the box size.')
        form.addParam('boxSize', IntParam,
                      label='Box size (px)',
                      validator=GE(0),
                      allowsNull=False,
                      help='Box size, in pixels, of the reconstruction. Note that this is independent of the '
                           'box size used to refine the particle. This allows the user to construct a 3D map of '
                           'arbitrary size to gain an overview of the structure surrounding the particle. A '
                           'sufficiently large box size also allows more of the high-frequency signal to be '
                           'captured that has been delocalized by the CTF.')
        form.addParam('croppedBoxSize', IntParam,
                      label="Cropped box size (px)",
                      allowsNull=True,
                      help='The resulting pseudo subtomograms are cropped to this size. A smaller box size '
                           'allows the (generally expensive) refinement using relion_refine to proceed more rapidly.')

    # -------------------------- STEPS functions ------------------------------
    def convertInputStep(self):
        self.genInStarFile(are2dParticles=self.getInputParticles().are2dStacks())

    # --------------------------- UTILS functions -----------------------------
    def _genCommonExtractAndRecCmd(self):
        cmd = [f'--p {self._getExtraPath(IN_PARTICLES_STAR)}',
               f'--t {self._getExtraPath(IN_TOMOS_STAR)}',
               f'--o {self._getExtraPath()}', f"--b {self.boxSize.get()}",
               f"--crop {self.croppedBoxSize.get()}",
               f"--bin {self.binningFactor.get():.1f}",
               f"--j {self.numberOfThreads.get()}",
               self._genExtraParamsCmd()]
        inParticles = self.getInputParticles()
        if type(inParticles) is RelionSetOfPseudoSubtomograms:
            trajectoriesFile = inParticles.getTrajectoriesStar()
            if trajectoriesFile:
                cmd.append(f'--mot {trajectoriesFile}')
        return ' '.join(cmd)
