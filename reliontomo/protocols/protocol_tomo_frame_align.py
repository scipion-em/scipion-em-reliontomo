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
from pyworkflow.protocol import IntParam, BooleanParam, GE, LE, FloatParam, EnumParam
from reliontomo import Plugin
from reliontomo.protocols.protocol_base_per_part_per_tilt import ProtRelionPerParticlePerTiltBase
from reliontomo.protocols.protocol_base_relion import IS_RELION_50
from reliontomo.utils import getProgram


class alignModels(Enum):
    tiltImages = 0
    onlyParticles = 1


class deformationModels(Enum):
    linear = 0
    spline = 1
    fourier = 2


class ProtRelionTomoFrameAlign(ProtRelionPerParticlePerTiltBase):
    """Relion has also implemented the analogous to Bayesian polishing in 2D for tomography.
     This procedure refines the projections that map 3D space onto the images of the
     tilt series. Optionally, the beam-induced motion trajectories of the particles
     and deformations can also be estimated.\n

     Each projection is optimised using the full 5 degrees of freedom: the assumption of
     a common tilt axis is abandoned. Even if no beam-induced motion is estimated,
     the (in that case static) 3D particle-positions are also optimised by the program.
     This is because those 3D positions cannot be assumed to be known in advance, since
     they have only been inferred from the observed 2D particle-positions in the
     individual tilt images. Therefore, this program always looks for the optimal set
     of 3D positions and projections.
    """

    _label = 'Bayesian polishing' if IS_RELION_50 else 'Frame alignment'

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        super()._defineParams(form)
        form.addSection(label='Polish')
        super()._insertBoxSizeForEstimationParam(form)
        form.addParam('maxPosErr', IntParam,
                      label='Max position error (px)',
                      default=5,
                      allowsNull=False,
                      validators=[GE(0), LE(64)],
                      help="Maximal assumed error in the initial 2D particle-positions (distances between the "
                           "projected 3D positions and their true positions in the images), given in pixels.")
        form.addParam('alignByShift', BooleanParam,
                      default=False,
                      label='Align by shift only?',
                      help='If set to Yes, tilt series projection shifts are refined based on cross-correlation. '
                           'Useful for very badly aligned frames. No iterative optimisation.')
        form.addParam('alignmentModel', EnumParam,
                      display=EnumParam.DISPLAY_HLIST,
                      choices=[choice.name for choice in alignModels],
                      label='Alignment model',
                      condition='alignByShift',
                      default=alignModels.tiltImages.value,
                      help='If set to "Only particles", it estimates rigid shift by aligning only the particles '
                           'instead of by predicting entire micrographs. In this case, only misalignments smaller than '
                           'half the box size of the particle can be corrected.')

        form.addSection(label='Motion & deformations')
        form.addParam('fitPerParticleMotion', BooleanParam,
                      default=False,
                      label='Fit per particle motion?',
                      help='If set to Yes, then the subtomogram version of Bayesian polishing will be used to fit '
                           'per-particle (3D) motion tracks, besides the rigid part of the motion in the tilt series.')
        group = form.addGroup('Per particle motion', condition='fitPerParticleMotion')
        group.addParam('sigmaVel', FloatParam,
                       label='Sigma for velocity (Å/dose)',
                       default=0.2,
                       validators=[GE(0.1), LE(10)],
                       help='The expected amount of motion (i.e. the std. deviation of particle positions in Angstroms '
                            'after 1 electron per A^2 of radiation).')
        group.addParam('sigmaDiv', IntParam,
                       label='Sigma for divergence (Å)',
                       default=5000,
                       validators=[GE(0), LE(10000)],
                       help='The expected spatial smoothness of the particle trajectories in angstroms (a greater '
                            'value means spatially smoother motion.')
        group.addParam('doGaussianDecay', BooleanParam,
                       label="Use Gaussian decay?",
                       default=False,
                       help='If set to Yes, then it is assumed that the correlation of the velocities of two particles '
                            'decays as a Gaussian over their distance, instead of as an exponential. This will produce '
                            'spatially smoother motion and result in a shorter program runtime.')
        if not IS_RELION_50:
            form.addParam('estimate2dDeformations', BooleanParam,
                          label='Estimate 2D deformations?',
                          condition='fitPerParticleMotion',
                          default=False,
                          help='If set to Yes, then the subtomogram version of Bayesian polishing will be used to fit '
                               'per-particle (3D) motion tracks, besides the rigid part of the motion in the tilt series.')
            group = form.addGroup('2D deformation estimation', condition='estimate2dDeformations')
            group.addParam('nHorizSamplingPts', IntParam,
                           label='Horizontal sampling points',
                           default=3,
                           validators=[GE(0), LE(10)],
                           help='Number of horizontal sampling points for the deformation grid.')
            group.addParam('nVertSamplingPts', IntParam,
                           label='Vertical sampling points',
                           default=3,
                           validators=[GE(0), LE(10)],
                           help='Number of vertical sampling points for the deformation grid.')
            group.addParam('deformationModel', EnumParam,
                           choices=[choice.name for choice in deformationModels],
                           label='Alignment model',
                           default=deformationModels.spline.value,
                           help='Type of model to use (linear, spline or Fourier).')
            group.addParam('deformationRegularisation', FloatParam,
                           label='Deformation regularisation scale',
                           default=0,
                           validators=[GE(0), LE(1)],
                           help='This is the strength of the deformation regularizer')
            group.addParam('refineDefPerFrame', BooleanParam,
                           label='Refine deformations per frame?',
                           default=False,
                           help='If set to Yes, it models deformations per tilt frame instead of per tilt series.')
        self._defineExtraParams(form)
        form.addParallelSection(threads=4, mpi=1)

    # -------------------------- INSERT steps functions -----------------------
    def _insertAllSteps(self):
        self._insertFunctionStep(self.convertInputStep)
        self._insertFunctionStep(self._relionTomoFrameAlign)
        self._insertFunctionStep(self.createOutputStep)

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
        # Polish
        cmd += '--r %i ' % self.maxPosErr.get()
        if self.alignByShift.get():
            cmd += '--shift_only '
            if self.alignmentModel.get() == alignModels.onlyParticles.value:
                cmd += '--shift_only_by_particles '

        # Motion
        if self.fitPerParticleMotion.get():
            cmd += '--motion '
            cmd += '--s_vel %.1f ' % self.sigmaVel.get()
            cmd += '--s_div %i ' % self.sigmaDiv.get()
            if self.doGaussianDecay.get():
                cmd += '--sq_exp_ker '
            if not IS_RELION_50:
                if self.estimate2dDeformations.get():
                    cmd += '--deformation '
                    cmd += '--def_w %i ' % self.nHorizSamplingPts.get()
                    cmd += '--def_h %i ' % self.nVertSamplingPts.get()
                    cmd += '--def_model %i ' % self.deformationModel.get()
                    cmd += '--def_reg %.2f ' % self.deformationRegularisation.get()
                    if self.refineDefPerFrame.get():
                        cmd += '--per_frame_deformation '

        return cmd
