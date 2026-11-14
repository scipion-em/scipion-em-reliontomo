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
    """
    Refines tomographic tilt-series projections and particle positions using
    Bayesian optimization strategies adapted for cryo-electron tomography.
    The protocol improves the consistency between the reconstructed 3D
    particle coordinates and their corresponding 2D projections across the
    tilt series, enabling more accurate subtomogram reconstruction and
    downstream structural interpretation.

    AI Generated:

    Bayesian Polishing and Frame Alignment (ProtRelionTomoFrameAlign) —
    User Manual

        Overview

        This protocol performs refinement of tomographic projection geometry
        and particle positioning in cryo-electron tomography datasets using
        Relion-based Bayesian optimization procedures. Its main objective is
        to improve the agreement between experimental tilt images and the
        reconstructed three-dimensional particle model by refining how each
        projection maps 3D space onto the observed images.

        In practical biological workflows, this refinement step is commonly
        used after initial subtomogram alignment and reconstruction. The
        protocol is particularly valuable when residual alignment errors,
        beam-induced motion, or local deformations limit the achievable
        resolution of subtomogram averages. By improving projection accuracy,
        the protocol can significantly enhance structural detail and overall
        map quality.

        Unlike simpler rigid alignment procedures, this approach optimizes
        the full geometric relationship between particles and tilt images.
        The refinement does not assume that particle positions are perfectly
        known beforehand. Instead, particle coordinates and projection
        parameters are refined simultaneously to produce a more physically
        consistent tomographic model.

        Projection Refinement and Position Optimization

        The protocol refines the orientation and positioning of each tilt
        image independently, allowing corrections that go beyond the
        assumption of a common tilt axis. This flexibility is important in
        real experimental datasets where small deviations in microscope
        geometry, stage movement, or alignment inaccuracies accumulate across
        the tilt series.

        From a biological perspective, accurate projection refinement is
        essential because small geometric inconsistencies can blur structural
        features during subtomogram averaging. Even modest improvements in
        projection accuracy may lead to clearer densities, improved secondary
        structure visibility, and better interpretation of flexible or
        heterogeneous complexes.

        The protocol also optimizes the three-dimensional particle positions.
        This is biologically important because particle coordinates inferred
        from noisy 2D observations are often imperfect. Refining these
        coordinates improves the consistency of the reconstructed tomogram
        and reduces errors propagated into downstream averaging workflows.

        Alignment by Shift Only

        An optional simplified alignment strategy allows refinement based
        exclusively on translational shifts. This mode is primarily intended
        for severely misaligned datasets where a rapid correction is needed
        before more advanced optimization can succeed.

        In biological practice, this option may be useful for problematic
        tilt series affected by poor frame alignment, stage instability, or
        preprocessing artifacts. Since the method relies on cross-correlation
        rather than iterative geometric refinement, it is generally faster
        but less flexible than the full optimization strategy.

        The protocol also allows particle-based shift estimation, where
        alignment focuses primarily on particle regions instead of entire
        micrographs. This can improve robustness when background signal,
        contamination, or large non-particle features interfere with global
        alignment. However, the method is best suited for relatively small
        misalignments and may fail when shifts exceed the particle box size.

        Per-Particle Motion Estimation

        One of the most biologically important features of this protocol is
        the ability to estimate per-particle motion trajectories. During
        electron irradiation, particles embedded in vitreous ice can move in
        complex ways that vary locally across the specimen. These motions are
        particularly relevant in tomography because the cumulative electron
        dose is distributed across many tilted projections.

        The protocol models these trajectories in three dimensions using
        Bayesian regularization principles. The resulting corrections improve
        coherence between tilt images and help recover high-resolution
        structural information that would otherwise be blurred by motion.

        The expected amount of particle motion can be controlled through
        velocity-related regularization parameters. Smaller values impose
        tighter constraints and assume relatively stable particles, whereas
        larger values permit stronger motion variability. Selecting
        biologically reasonable values is important because excessive freedom
        may lead to overfitting noise rather than true motion.

        Spatial smoothness constraints regulate how neighboring particles are
        expected to move relative to each other. In biological samples,
        nearby particles embedded within the same ice region often exhibit
        partially correlated motion due to local mechanical deformation.
        Increasing spatial smoothness produces more coherent motion fields,
        while lower smoothness allows more localized particle behavior.

        Gaussian or exponential spatial decay models can also be selected to
        describe how motion correlations decrease with distance. Gaussian
        decay generally produces smoother motion fields and may reduce
        computational complexity in large datasets.

        Estimation of Local Deformations

        In some workflow configurations, the protocol can additionally model
        two-dimensional local deformations across tilt images. These
        deformations represent non-rigid distortions caused by beam-induced
        specimen bending, local charging, or stage-related instabilities.

        From a biological perspective, deformation correction becomes
        especially important in thick specimens, crowded cellular
        environments, or datasets acquired under challenging imaging
        conditions. Correcting these effects can improve the consistency of
        particle alignment across the tilt series and increase the quality of
        subtomogram averages.

        The deformation field is sampled across the image using a configurable
        grid. Denser sampling allows finer local corrections but increases
        computational complexity and the risk of overfitting. In many
        biological applications, moderate sampling densities provide a good
        balance between flexibility and robustness.

        Different deformation models may be available, including linear,
        spline, or Fourier-based representations. Simpler models are often
        sufficient for mild distortions, whereas spline or Fourier approaches
        can better capture complex spatial variations in larger specimens.

        A regularization parameter controls the smoothness and stability of
        the deformation field. Stronger regularization suppresses unrealistic
        local distortions, while weaker regularization allows more flexible
        corrections. Excessively weak regularization may introduce artifacts
        that lack biological meaning.

        Per-frame deformation refinement can also be enabled, allowing the
        protocol to estimate independent deformation behavior for each tilt
        frame rather than a single deformation model for the entire tilt
        series. This option is computationally demanding but may improve
        results in highly dynamic or unstable acquisitions.

        Inputs and Workflow Considerations

        The protocol typically operates on aligned tilt-series data together
        with particle coordinates and subtomogram metadata derived from
        previous tomography processing steps. Reliable initial geometry is
        important because the refinement assumes that the dataset is already
        approximately consistent.

        In practice, users often apply this protocol after obtaining an
        initial subtomogram reconstruction and before performing final
        high-resolution averaging or classification. The refinement may be
        iterated multiple times as part of progressive optimization
        workflows.

        Biological users should ensure that tilt-series preprocessing,
        contrast transfer correction, and particle extraction have been
        performed carefully before running this protocol. Severe artifacts,
        missing tilt information, or inaccurate initial coordinates may limit
        the effectiveness of the refinement.

        Outputs and Biological Interpretation

        The protocol produces refined geometric parameters, corrected
        particle trajectories, and improved particle positions for the tilt
        series. These outputs are subsequently used in subtomogram averaging
        and refinement workflows to generate higher-quality reconstructions.

        Improved geometric consistency often translates directly into better
        map interpretability. Biological structures may display sharper
        secondary structure features, improved membrane definition, or more
        reliable density continuity after refinement.

        Motion estimation results can also provide indirect information about
        specimen stability and imaging quality. Large or irregular motion
        patterns may indicate problematic regions of the sample, excessive
        beam-induced deformation, or challenging acquisition conditions.

        Practical Recommendations

        For most datasets, it is advisable to begin with conservative motion
        estimation parameters and moderate regularization values. Excessively
        flexible refinement may fit noise instead of meaningful motion or
        deformation.

        When processing relatively stable samples with good initial
        alignments, standard projection refinement without deformation
        modeling is often sufficient. More advanced motion and deformation
        estimation should typically be reserved for challenging datasets,
        thick specimens, or high-resolution optimization workflows.

        Users should visually inspect reconstruction quality before and after
        refinement to confirm that the corrections improve biological
        interpretability rather than introducing artifacts. Careful parameter
        tuning and iterative refinement are often necessary to achieve the
        best results in demanding cryo-electron tomography studies.
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
        form.addParallelSection(threads=0, mpi=1)

    # -------------------------- INSERT steps functions -----------------------
    def _insertAllSteps(self):
        self._insertFunctionStep(self.convertInputStep, needsGPU=False)
        self._insertFunctionStep(self._relionTomoFrameAlign, needsGPU=False)
        self._insertFunctionStep(self.createOutputStep, needsGPU=False)

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
