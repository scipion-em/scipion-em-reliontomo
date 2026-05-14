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
from os import remove
from os.path import abspath, exists
from pyworkflow.object import Boolean
from pyworkflow.protocol import LEVEL_ADVANCED, IntParam, StringParam, BooleanParam, \
    EnumParam, PathParam, FloatParam, LEVEL_NORMAL, GE, LE
from reliontomo.constants import ANGULAR_SAMPLING_LIST
from reliontomo.protocols.protocol_base_relion import ProtRelionTomoBase, IS_RELION_50


class ProtRelionRefineBase(ProtRelionTomoBase):
    """
    Provides the common refinement framework used for subtomogram
    reconstruction, initial model generation, and iterative auto-refinement
    workflows within RELION tomography environments.

    AI Generated:

    RELION Tomography Refinement Base (ProtRelionRefineBase) — User Manual
        Overview

        This protocol serves as the central refinement foundation for
        tomography-based reconstruction workflows in RELION. Its purpose is
        to prepare, control, and optimize iterative refinement procedures for
        pseudo-subtomograms and related particle representations generated
        from cryo-electron tomography experiments. The protocol establishes a
        consistent environment for both ab initio model generation and
        high-resolution auto-refinement while supporting a broad range of
        experimental conditions and computational infrastructures.

        In biological cryo-ET studies, refinement is the stage where noisy
        particle observations are progressively aligned and combined into a
        coherent three-dimensional reconstruction. This process is essential
        for improving structural interpretability, separating heterogeneous
        conformations, and obtaining biologically meaningful density maps.
        The protocol is therefore commonly used after particle extraction and
        before downstream interpretation steps such as classification, atomic
        modeling, or focused analysis.

        Input Data and Reconstruction Context

        The protocol operates on pseudo-subtomograms generated from tilt-series
        data. These inputs represent intermediate reconstruction objects that
        approximate the original two-dimensional information while remaining
        compatible with the RELION refinement framework. They are designed to
        preserve alignment and contrast-transfer information required for
        iterative Bayesian optimization.

        From a biological perspective, the quality of the input particles
        strongly influences refinement stability. Accurate particle picking,
        reliable tomogram alignment, and appropriate CTF estimation are all
        critical prerequisites. Poorly aligned subtomograms or highly
        heterogeneous particle populations may lead to unstable refinements or
        loss of structural detail.

        The protocol also supports workflows involving both two-dimensional
        particle stacks and fully three-dimensional subtomographic data,
        allowing flexibility across different tomography processing strategies.

        CTF Correction and Optical Considerations

        Contrast transfer function correction plays a central role in cryo-EM
        and cryo-ET refinement because it compensates for microscope-induced
        phase distortions. This protocol allows refinement workflows to apply
        CTF correction during projection matching and reconstruction, improving
        the recovery of high-resolution information.

        In some experimental situations, especially when low-resolution CTF
        behavior is unreliable, users may choose to suppress correction near
        the first CTF peak. This can stabilize refinements for difficult data,
        although it may also reduce sensitivity to low-frequency information.
        For most well-behaved datasets, standard CTF correction remains the
        preferred strategy.

        Refinement and Regularisation Strategy

        The protocol supports iterative Bayesian refinement approaches in
        which experimental evidence is balanced against prior expectations.
        Regularisation controls this balance and strongly affects refinement
        behavior. Lower regularisation values tend to enforce smoother and
        more conservative reconstructions, while larger values increase the
        influence of experimental observations and may reveal additional
        structural detail.

        In practical biological applications, moderate regularisation values
        are often preferred during early refinement stages to avoid overfitting
        noise. As the refinement stabilizes and particle quality improves,
        stronger weighting of the experimental data may become beneficial.

        The framework also supports stochastic gradient optimization strategies
        suitable for initial model generation and heterogeneous datasets. These
        approaches are particularly valuable when the initial structure is
        unknown or when multiple structural states coexist within the dataset.

        Classification and Structural Heterogeneity

        The protocol can support workflows involving multiple structural
        classes. In biological cryo-ET studies, heterogeneous conformations
        are extremely common due to molecular flexibility, compositional
        variability, or differences in biochemical state.

        Multi-class refinement can separate these structural populations and
        prevent incompatible particles from degrading reconstruction quality.
        In practice, users often employ multiple classes to identify damaged
        particles, contaminants, assembly intermediates, or distinct functional
        conformations. Homogeneous datasets may refine successfully with a
        single class, whereas highly variable systems typically benefit from
        more extensive classification strategies.

        Symmetry and Structural Constraints

        Symmetry handling is one of the most biologically important aspects of
        refinement because it directly affects the interpretation of the final
        reconstruction. Correct symmetry application can dramatically improve
        signal quality by averaging equivalent structural views. However,
        imposing incorrect symmetry may artificially distort biologically
        meaningful asymmetry.

        For uncertain systems, it is often advisable to begin refinement in
        asymmetric mode and only apply symmetry after validating the particle
        architecture. This approach is especially important for flexible
        assemblies, partially occupied complexes, and macromolecules with
        pseudo-symmetry.

        The protocol also supports workflows in which symmetry determination
        is delayed until later refinement stages. This strategy can improve
        robustness during early optimization and reduce the risk of convergence
        toward incorrect symmetric solutions.

        Particle Masking and Solvent Treatment

        Proper masking is critical for stable refinement because it determines
        which regions contribute most strongly to alignment and reconstruction.
        Circular particle masks are commonly used to isolate the molecular
        signal while suppressing solvent noise.

        For compact particles, standard spherical masking is usually sufficient.
        More complex systems, including elongated assemblies or membrane
        proteins with flexible regions, may require careful adjustment of mask
        diameter to avoid excluding biologically relevant density.

        Solvent flattening and non-negative density enforcement further improve
        reconstruction stability by reducing unrealistic density fluctuations.
        These constraints are particularly useful during early iterations or
        when processing noisy subtomograms.

        Angular Sampling and Alignment Precision

        Angular and translational sampling parameters determine how broadly the
        protocol explores possible particle orientations and positions. Wider
        searches improve robustness for poorly aligned datasets but increase
        computational cost. Narrower searches improve efficiency and precision
        once approximate alignments are already available.

        Prior angular information can also be incorporated into refinement
        workflows. This is especially valuable in tomography datasets where
        particles may exhibit biologically constrained orientations, such as
        membrane-associated complexes, filaments, or surface-bound assemblies.

        The careful use of priors can significantly stabilize refinement while
        preserving meaningful structural variability.

        GPU Acceleration and Computational Scaling

        Tomographic refinement is computationally demanding because it involves
        repeated alignment and reconstruction of large particle datasets. The
        protocol therefore supports parallel execution strategies, GPU
        acceleration, memory optimization, and advanced disk-management modes.

        GPU acceleration substantially reduces refinement time and is strongly
        recommended for large-scale datasets or high-resolution studies. Memory
        management options allow users to optimize performance depending on
        available hardware resources, storage speed, and network configuration.

        On high-performance computing systems, balancing MPI processes, thread
        usage, GPU distribution, and memory consumption is essential for stable
        and efficient execution.

        Output Interpretation

        The protocol produces refined particle metadata, reconstruction maps,
        and iteration-dependent refinement information that can be used for
        downstream processing. These outputs form the basis for subsequent
        classification, post-processing, local refinement, focused masking,
        and biological interpretation.

        During iterative refinement, users should monitor convergence behavior,
        map quality, and consistency between iterations. Sudden changes in map
        appearance, unstable angular assignments, or unexpected resolution
        behavior may indicate overfitting, heterogeneity, or parameter
        misconfiguration.

        Practical Recommendations

        For most biological cryo-ET projects, it is advisable to begin with
        conservative refinement settings and gradually increase refinement
        complexity as stability improves. Initial reconstructions often benefit
        from broad angular searches, moderate regularisation, and careful
        masking strategies.

        GPU acceleration should be enabled whenever compatible hardware is
        available. Users should also verify that reference maps, masks, and
        particles share compatible dimensions and sampling rates before
        starting refinement.

        In heterogeneous datasets, introducing multiple classes early in the
        workflow often improves refinement robustness and prevents poor-quality
        particles from dominating the reconstruction. Symmetry should only be
        applied when biologically justified.

        Final Perspective

        Tomographic refinement is not simply a numerical optimization process
        but a biologically meaningful reconstruction strategy that transforms
        noisy particle observations into interpretable structural information.
        Successful refinement depends on balancing computational precision,
        prior biological knowledge, experimental quality, and careful parameter
        selection. When used appropriately, this protocol provides a powerful
        foundation for high-resolution cryo-electron tomography analysis.
    """

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.alignmentAsPriors = Boolean(False)

    # -------------------------- DEFINE param functions -----------------------

    # I/O PARAMS -------------------------------------------------------------------------------------------------------
    def _defineIOParams(self, form):
        self._defineCommonInputParams(form)
        self._insertBinThreadsParam(form)

    # CTF PARAMS -------------------------------------------------------------------------------------------------------
    @staticmethod
    def _defineCTFParams(form):
        form.addSection(label='CTF')
        form.addParam('doCTF', BooleanParam,
                      default=True,
                      label='Do CTF-correction?',
                      help='If set to Yes, CTFs will be applied to the projections of the map. '
                           'Note that input particles should contains CTF parameters.')
        form.addParam('ignoreCTFUntilFirstPeak', BooleanParam,
                      default=False,
                      label='Ignore CTFs until first peak?',
                      help='If set to Yes, then CTF-amplitude correction will only be performed from the first peak '
                           'of each CTF onward. This can be useful if the CTF model is inadequate at the lowest '
                           'resolution. Still, in general using higher amplitude contrast on the CTFs (e.g. 10-20%) '
                           'often yields better results. Therefore, this option is not generally recommended.')

    # OPTIMIZATION PARAMS ----------------------------------------------------------------------------------------------
    @staticmethod
    def _insertOptimisationSection(form):
        form.addSection(label='Optimisation')

    @staticmethod
    def _insertVdamMiniBatchesParam(form):
        form.addParam('nVdamMiniBatches', IntParam,
                      allowsNull=False,
                      default=200,
                      label='Number of VDAM mini-batches',
                      help="How many iterations (i.e. mini-batches) to perform with the VDAM ((variable metric "
                           "gradient descent with adaptive moments) algorithm?. Using 200 (default) has given good "
                           "results for many data sets. Using 100 will run faster, at the expense of some quality in "
                           "the results.")

    @staticmethod
    def _insertRegularisationParam(form, isCl3d=False):
        form.addParam('regularisation', FloatParam,
                      default=1 if isCl3d else 4,
                      validators=[GE(0)],
                      label='Regularisation parameter T',
                      help="Bayes law strictly determines the relative weight between the contribution of the "
                           "experimental data and the prior. However, in practice one may need to adjust this weight "
                           "to put slightly more weight on the experimental data to allow optimal results. If it's set "
                           "to 0, no regularisation will be applied. Values greater than 1 for this regularisation "
                           "parameter (T in the JMB2011 paper) put more weight on the experimental data. Values "
                           "around 2-4 have been observed to be useful for 3D initial model calculations.")

    @staticmethod
    def _insertNumOfClassesParam(form):
        form.addParam('numberOfClasses', IntParam,
                      default=1,
                      validators=[GE(1)],
                      label='Number of classes',
                      help='The number of classes (K) for a multi-reference ab initio SGD refinement. These classes '
                           'will be made in an unsupervised manner, starting from a single reference in the initial '
                           'iterations of the SGD, and the references will become increasingly dissimilar during the '
                           'inbetween iterations. Sometimes, using more than one class may help in providing a ‘sink’ '
                           'for sub-optimal particles that may still exist in the data set. In this case, which is '
                           'quite homogeneous, a single class should work just fine.')

    @staticmethod
    def _insertBlushRegParam(form):
        form.addParam('doBlushReg', BooleanParam,
                      default=False,
                      label='Use blush regularisation',
                      help='If set to Yes, relion_refine will use a neural network to perform regularisation by '
                           'denoising at every iteration, instead of the standard smoothness regularisation')

    @staticmethod
    def _insertMaskDiameterParam(form):
        form.addParam('maskDiameter', IntParam,
                      allowsNull=False,
                      validators=[GE(0)],
                      label='Circular mask diameter (Å)',
                      help='The experimental images will be masked with a soft circular mask with this diameter. '
                           'Make sure this radius is not set too small because that may mask away part of the signal! '
                           'If set to a value larger than the image size no masking will be performed. \n '
                           'The same diameter will also be used for a spherical mask of the reference structures if '
                           ' no user-provided mask is specified.')

    @staticmethod
    def _insertZeroMaskParam(form):
        form.addParam('zeroMask', BooleanParam,
                      label='Mask particles with zeros?',
                      default=True,
                      help="If set to Yes, then in the individual particles, the area outside a circle with the "
                           "radius of the particle will be set to zeros prior to taking the Fourier transform.\n\nThis "
                           "will remove noise and therefore increase sensitivity in the alignment and classification.\n"
                           "\nHowever, it will also introduce correlations between the Fourier components that are not "
                           "modelled. When set to No, then the solvent area is filled with random noise, which "
                           "prevents introducing correlations.\n\nHigh-resolution refinements (e.g. ribosomes or other "
                           "large complexes in 3D auto-refine) tend to work better when filling the solvent area with "
                           "random noise (i.e. setting this option to No), refinements of smaller complexes and most "
                           "classifications go better when using zeros (i.e. setting this option to Yes).")

    @staticmethod
    def _insertFlattenSolventParam(form):
        form.addParam('flattenSolvent', BooleanParam,
                      default=True,
                      label='Flatten and enforce non-negative solvent?',
                      help="If set to Yes, the job will apply a spherical mask and enforce all values in the "
                           "reference to be non-negative.")

    @staticmethod
    def _insertSymmetryParam(form, helpMsg):
        form.addParam('symmetry', StringParam,
                      label='Symmetry group',
                      default='C1',
                      help=helpMsg)

    @staticmethod
    def _insertDoInC1AndApplySymLaterParam(form):
        form.addParam('doInC1AndApplySymLater', BooleanParam,
                      default=True,
                      label='Run in C1 and apply symmetry later?',
                      help="If set to Yes, the gradient-driven optimisation is run in C1 and the symmetry orientation "
                           "is searched and applied later. If set to No, the entire optimisation is run in the "
                           "symmetry point group indicated above.")

    @staticmethod
    def _insertPriorWidthParam(form):
        form.addParam('priorWidthTiltAngle', IntParam,
                      default=-1,
                      label='Prior width on tilt angle (deg)',
                      help='The width of the prior on the tilt angle: angular searches will be +/-3 times this value. '
                           'Tilt priors will be defined when particles have been picked as filaments, on spheres or '
                           'on manifolds. Setting this width to a negative value will lead to no prior being used on '
                           'the tilt angle.')

    # COMPUTE PARAMS ---------------------------------------------------------------------------------------------------
    def _defineComputeParams(self, form, isOnlyClassif=False):
        form.addSection(label='Compute')
        form.addParam('parallelDiscIO', BooleanParam,
                      default=True,
                      label='Use parallel disc I/O?',
                      help="If set to Yes, all MPI followers will read their own images from disc. Otherwise, only "
                           "the leader will read images and send them through the network to the followers. Parallel "
                           "file systems like gluster of fhgfs are good at parallel disc I/O. NFS may break with "
                           "many followers reading in parallel. If your datasets contain particles with different "
                           "box sizes, you have to say Yes.")
        form.addParam('pooledSubtomos', IntParam,
                      default=1,
                      validators=[GE(1)],
                      label='Number of pooled particles',
                      help="Particles are processed in individual batches by MPI followers. During each batch, a "
                           "stack of particle images is only opened and closed once to improve disk access times.\n\n"
                           "All particle images of a single batch are read into memory together. The size of these "
                           "batches is at least one particle per thread used. This parameter controls how many "
                           "particles are read together for each thread. If it is set to 3 and one uses 8 threads, "
                           "batches of 3x8=24 particles will be read together.\n\nThis may improve performance on "
                           "systems where disk access, and particularly metadata handling of disk access, is a "
                           "problem. It has a modest cost of increased RAM usage.")
        if not IS_RELION_50:
            form.addParam('skipGridding', BooleanParam,
                          default=True,
                          label='Skip griding?',
                          help='Skip gridding in the Maximization step in the Expectation-Maximization algorithm. '
                               'If this option is set to Yes, more memory will be consumed during the protocol '
                               'execution, but it will be faster.')
        form.addParam('allParticlesRam', BooleanParam,
                      default=False,
                      label='Pre-read all particles into RAM?',
                      help='If set to Yes, all particle images will be read into computer memory, which will greatly '
                           'speed up calculations on systems with slow disk access. However, one should of course '
                           'be careful with the amount of RAM available. Because particles are read in float-precision,'
                           ' it will take ( N * box_size * box_size * 4 / (1024 * 1024 * 1024) ) Giga-bytes to'
                           ' read N particles into RAM. For 100 thousand 200x200 images, that becomes 15Gb, or 60 Gb'
                           ' for the same number of 400x400 particles. Remember that running a single MPI follower '
                           'on each node that runs as many threads as available cores will have access to all '
                           'available RAM.\n'
                           'If parallel disc I/O is set to No, then only the leader reads all particles into RAM and '
                           'sends those particles through the network to the MPI followers during the refinement '
                           'iterations.')
        form.addParam('scratchDir', PathParam,
                      label='Copy particles to scratch directory',
                      help='If a directory is provided here, then the job will create a sub-directory in it called'
                           ' relion_volatile. If that relion_volatile directory already exists, it will be wiped. '
                           'Then, the program will copy all input particles into a large stack inside the '
                           'relion_volatile '
                           'subdirectory. Provided this directory is on a fast local drive (e.g. an SSD drive), '
                           'processing in all the iterations will be faster. If the job finishes correctly, '
                           'the relion_volatile directory will be wiped. If the job crashes, you may want to remove it '
                           'yourself.')
        form.addParam('combineItersDisc', BooleanParam,
                      default=False,
                      label='Combine iterations through disc?',
                      help='If set to Yes, at the end of every iteration all MPI followers will write out a large '
                           'file with their accumulated results. The MPI leader will read in all these files, '
                           'combine them all, and write out a new file with the combined results. All MPI salves '
                           'will then read in the combined results. This reduces heavy load on the network, '
                           'but increases load on the disc I/O. This will affect the time it takes between the '
                           'progress-bar in the expectation step reaching its end (the mouse gets to the cheese) '
                           'and the start of the ensuing maximisation step. It will depend on your system setup'
                           ' which is most efficient.')

    # ADDITIONAL PARAMS ------------------------------------------------------------------------------------------------
    @staticmethod
    def _defineAdditionalParams(form):
        form.addSection(label='Additional')
        form.addParam('keepOnlyLastIterFiles',
                      BooleanParam,
                      default=False,
                      expertLevel=LEVEL_ADVANCED,
                      label="Keep only files from last iteration?",
                      help="If Yes is chosen, only the files which correspond to the last iteration will be saved "
                           "in the protocol's extra directory. Otherwise, files corresponding to each iteration "
                           "will be kept.")
        form.addParam('oversampling', IntParam,
                      default=1,
                      expertLevel=LEVEL_ADVANCED,
                      label="Over-sampling",
                      help="Adaptive oversampling order to speed-up calculations (0=no oversampling, 1=2x, 2=4x, etc)")

        # Generate priors fields
        if not IS_RELION_50:
            form.addParam('alignmentAsPriors', BooleanParam, default=False,
                          expertLevel=LEVEL_ADVANCED,
                          label='Consider alignment as priors?',
                          help='If set to Yes, then alignment information from '
                               'input particles will be considered as PRIORS. This '
                               'option can be used to do restricted local '
                               'search within a range centered around those priors.')

        ProtRelionTomoBase._defineExtraParams(form, addAdditionalSection=False)
        form.addParallelSection(threads=0, mpi=1)

    # ANGULAR SAMPLING PARAMS ------------------------------------------------------------------------------------------
    @staticmethod
    def _insertAngularCommonParams(form, expertLevel=LEVEL_NORMAL, angSampling=2, offsetRange=5,
                                   offsetStep=1, condition=True):
        form.addParam('angularSamplingDeg', EnumParam,
                      default=angSampling,
                      condition=condition,
                      choices=ANGULAR_SAMPLING_LIST,
                      label='Initial angular sampling interval (deg)',
                      expertLevel=expertLevel,
                      help='There are only a few discrete angular samplings possible because '
                           'we use the HealPix library to generate the sampling of the first '
                           'two Euler angles on the sphere. The samplings are approximate numbers '
                           'and vary slightly over the sphere.\n  Note that this will only be the '
                           'value for the first few iteration(s): the sampling rate will be '
                           'increased automatically after that.')
        form.addParam('offsetSearchRangePix', FloatParam,
                      default=offsetRange,
                      condition=condition,
                      label='Initial offset range (pix.)',
                      expertLevel=expertLevel,
                      validators=[GE(0), LE(30)],
                      help='Probabilities will be calculated only for translations in a circle '
                           'with this radius (in pixels). The center of this circle changes at '
                           'every iteration and is placed at the optimal translation for each '
                           'image in the previous iteration.')
        form.addParam('offsetSearchStepPix', FloatParam,
                      default=offsetStep,
                      condition=condition,
                      label='Initial offset step (pix.)',
                      expertLevel=expertLevel,
                      validators=[GE(0.1), LE(5)],
                      help='Translations will be sampled with this step-size (in pixels). '
                           'Translational sampling is also done using the adaptive approach. '
                           'Therefore, if adaptive=1, the translations will first be evaluated'
                           'on a 2x coarser grid.')

    @staticmethod
    def _insertGpuParams(form):
        form.addParam('doGpu', BooleanParam,
                      default=False,
                      label='Use GPU acceleration?',
                      help='If set to Yes, it will use available gpu resources for some calculations.')
        form.addParam('gpusToUse', StringParam,
                      condition='doGpu',
                      default='0',
                      label='GPUs to use:',
                      help='It can be used to provide a list of which GPUs (e. g. "0:1:2:3") to use. MPI-processes are '
                           'separated by ":", threads by ",". For example: "0,0:1,1:0,0:1,1"')

    # -------------------------- INSERT steps functions -----------------------

    def convertInputStep(self):
        if IS_RELION_50:
            self.genInStarFile(are2dParticles=self.getInputParticles().are2dStacks())
        else:
            self.genInStarFile(withPriors=self.alignmentAsPriors)

    # -------------------------- INFO functions -------------------------------
    def _validate(self):
        def _getDims(obj, is2d):
            """ Helper to get (x, y) or (x, y, z) based on dimensionality. """
            x, y, z = obj.getDimensions()
            return (x, y) if is2d else (x, y, z)

        errorMsg = super()._validate()
        sRateTol = 1e-2
        inParticles = self.getInputParticles()
        inParticlesSRate = inParticles.getSamplingRate()
        is2D = inParticles.are2dStacks()
        partDims = _getDims(inParticles, is2D)

        # 1. Check Reference Volume
        if not self.doResizeRef.get():
            refVolume = self.referenceVolume.get()
            if abs(inParticlesSRate - refVolume.getSamplingRate()) >= sRateTol:
                errorMsg.append(f'The introduced particles and the reference volume must have the same sampling rate:\n'
                                f'{inParticlesSRate:.3f} != {refVolume.getSamplingRate():.3f} Å/px. Consider resizing.')

            if _getDims(refVolume, is2D) != partDims:
                errorMsg.append(f'The dimensions of the reference volume and the particles must be the same.')

        # 2. Check Solvent Mask
        solventMask = self.solventMask.get()
        if solventMask:
            maskSRate = solventMask.getSamplingRate()
            if abs(inParticlesSRate - maskSRate) >= sRateTol:
                errorMsg.append(f'The introduced particles and the solvent mask must have the same sampling rate:\n'
                                f'{inParticlesSRate:.3f} != {maskSRate:.3f} Å/px')

            if _getDims(solventMask, is2D) != partDims:
                errorMsg.append(f'The dimensions of the solvent mask and the particles must be the same.')

        # 3. Resource checks (MPI & GPUs)
        if self.numberOfMpi.get() < 3:
            errorMsg.append('The number of MPIs must be at least 3.')

        gpus = self.gpusToUse.get().replace(' ', '')
        if len(gpus) > 1 and ':' not in gpus:
            errorMsg.append("Bad format detected in GPUs -> they should be separated by ':'.")

        return errorMsg

    def _warnings(self):
        pass


    # --------------------------- UTILS functions -----------------------------
    def _genBaseCommand(self, useOptimizationSet=False):
        cmd = ''
        cmd += self._genIOBaseCmd(useOptimizationSet=useOptimizationSet)  # I/O args
        cmd += self._genCTFBaseCmd()  # CTF args
        cmd += self._genOptimisationBaseCmd()  # Optimisation args
        cmd += self._genComputeBaseCmd()  # Compute args
        cmd += self._genAdditionalBaseCmd()  # Additional args
        return cmd

    def _genIOBaseCmd(self, useOptimizationSet=False):
        inRelionParticles = self.getInputParticles()
        if useOptimizationSet:
            # Use optimization set file
            self.info("Using optimization_set: %s" % inRelionParticles.filesMaster)
            cmd = '--ios %s ' % inRelionParticles.filesMaster
        else:
            cmd = '--i %s ' % self.getOutStarFileName()
            cmd += '--tomograms %s ' % inRelionParticles.getTomogramsStar()
            trajectoriesFile = inRelionParticles.getTrajectoriesStar()
            if trajectoriesFile:
                cmd += '--trajectories %s ' % trajectoriesFile
        cmd += '--o %s ' % (self._getExtraPath() + '/')  # If not, Relion will concatenate it directly as a prefix
        cmd += '--j %i ' % self.binThreads.get()
        return cmd

    def _genCTFBaseCmd(self):
        cmd = ''
        if self.doCTF.get():
            cmd += '--ctf '
        if self.ignoreCTFUntilFirstPeak.get():
            cmd += '--ctf_intact_first_peak '
        return cmd

    def _genOptimisationBaseCmd(self):
        return '--particle_diameter %i ' % self.maskDiameter.get()

    def _genComputeBaseCmd(self, onlyCl3d=False):
        cmd = ''
        if not self.parallelDiscIO.get():
            cmd += '--no_parallel_disc_io '
        cmd += '--pool %i ' % self.pooledSubtomos.get()
        if self.allParticlesRam.get():
            cmd += '--preread_images '
        if not self.combineItersDisc.get():
            cmd += '--dont_combine_weights_via_disc '
        if self.scratchDir.get():
            cmd += '--scratch_dir %s ' % self.scratchDir.get()
        if self.doGpu.get() and not onlyCl3d:
            cmd += '--gpu "%s" ' % self.gpusToUse.get()
        return cmd

    def _genAdditionalBaseCmd(self):
        cmd = '--oversampling %i ' % self.oversampling.get()
        cmd += self._genExtraParamsCmd()
        return cmd

    def _applyKeepIterFilesUserSelection(self):
        if self.keepOnlyLastIterFiles.get():
            self._cleanUndesiredFiles()

    def _cleanUndesiredFiles(self):
        """Remove all files generated by relion_classify 3d excepting the ones which
        correspond to the last iteration. Example for iteration 25:
        relion_it025_class002.mrc
        relion_it025_class001.mrc
        relion_it025_model.star
        relion_it025_sampling.star
        relion_it025_optimiser.star
        relion_it025_data.star
        """
        itPref = 'relion_it'
        clPref = 'class'
        starExt = '.star'
        mrcExt = '.mrc'
        # Classify calculations related files
        calcFiles = ['data', 'model', 'optimiser', 'sampling']
        for i in range(self._lastIter()):
            for calcFile in calcFiles:
                fn = abspath(self._getExtraPath('{}{:03d}_{}{}'.format(
                    itPref, i, calcFile, starExt)))
                if exists(fn):
                    remove(fn)
            # Classes related files
            for itr in range(1, self.nClasses + 1):
                fn = abspath(self._getExtraPath('{}{:03d}_{}{:03d}{}'.format(
                    itPref, i, clPref, itr, mrcExt)))
                if exists(fn):
                    remove(fn)
