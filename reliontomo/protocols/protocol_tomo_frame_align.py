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
from pyworkflow.protocol import IntParam, BooleanParam, GE, LE, FloatParam
from reliontomo import Plugin
from reliontomo.protocols.protocol_base_per_part_per_tilt import ProtRelionPerParticlePerTiltBase
from reliontomo.utils import getProgram
from tomo.protocols import ProtTomoBase


class ProtRelionTomoFrameAlign(ProtRelionPerParticlePerTiltBase, ProtTomoBase):
    """Tomo frame align"""

    _label = 'Tomo frame align'

    def __init__(self, **args):
        ProtRelionPerParticlePerTiltBase.__init__(self, **args)

    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        ProtRelionPerParticlePerTiltBase._defineParams(self, form)
        form.addSection(label='Polish')
        ProtRelionPerParticlePerTiltBase._insertBoxSizeForEstimationParam(form)
        form.addParam('maxPosErr', IntParam,
                      label='Max position error (pix)',
                      default=5,
                      allowsNull=False,
                      validators=[GE(0), LE(64)],
                      help="Maximal assumed error in the initial 2D particle-positions (distances between the "
                           "projected 3D positions and their true positions in the images), given in pixels.")
        form.addParam('doFlexAlign', BooleanParam,
                      label="Allow flexible alignment?",
                      default=False,
                      help="If set to No, only an optimal rigid shift will be applied to each frame (no iterative "
                           "optimisation).")
        form.addParam('doGlobalRigidAlign', BooleanParam,
                      label="Do global rigid shift alignment?",
                      default=False,
                      condition='not doFlexAlign',
                      help="If set to Yes, it estimates the rigid shift by aligning only the particles instead of by "
                           "predicting the entire micrographs.")
        form.addParam('doPolish', BooleanParam,
                      label="Fit per-particle motion?",
                      default=False,
                      condition='doFlexAlign',
                      help="If set to Yes, then the subtomogram version of Bayesian polishing will be used to fit "
                           "per-particle (3D) motion tracks, besides the rigid part of the motion in the tilt series.")
        form.addParam('sigmaVel', FloatParam,
                      label='Sigma for velocity (Å/dose)',
                      condition='doPolish',
                      default=0.2,
                      validators=[GE(0.1), LE(10)],
                      help="The expected amount of motion (i.e. the std. deviation of particle positions in Angstroms "
                           "after 1 electron per A^2 of radiation).")
        form.addParam('sigmaDiv', IntParam,
                      label='Sigma for divergence (Å)',
                      default=5000,
                      condition='doPolish',
                      validators=[GE(0), LE(10000)],
                      help="The expected spatial smoothness of the particle trajectories in angstroms (a greater value "
                           "means spatially smoother motion.")
        form.addParam('doGaussianDecay', BooleanParam,
                      label="Use Gaussian decay?",
                      default=False,
                      condition='doPolish',
                      help="If set to Yes, then it's assumed that the correlation of the velocities of two particles "
                           "decays as a Gaussian over their distance, instead of as an exponential. This will produce "
                           "spatially smoother motion and result in a shorter program runtime.")

    # -------------------------- INSERT steps functions -----------------------
    def _insertAllSteps(self):
        self._insertFunctionStep(self._relionTomoFrameAlign)
        # self._insertFunctionStep(self.createOutputStep)

    # -------------------------- STEPS functions ------------------------------
    def _relionTomoFrameAlign(self):
        Plugin.runRelionTomo(self, getProgram('relion_tomo_align', self.numberOfMpi.get()),
                             self._genTomoFrameAlignCmd(), numberOfMpi=self.numberOfMpi.get())

    # -------------------------- INFO functions -------------------------------
    def _validate(self):
        pass

    # --------------------------- UTILS functions -----------------------------
    def _genTomoFrameAlignCmd(self):
        cmd = self._genIOCommand()
        cmd += '--r %i ' % self.maxPosErr.get()
        if self.doFlexAlign.get():
            if self.doPolish.get():
                cmd += '--motion '
                cmd += '--s_vel %.1f ' % self.sigmaVel.get()
                cmd += '--s_div %i ' + self.sigmaDiv.get()
                if self.doGaussianDecay.get():
                    cmd += '--sq_exp_ker '
        else:
            cmd += '--shift_only '
            if self.doGlobalRigidAlign.get():
                cmd += '--shift_only_by_particles '

        return cmd
