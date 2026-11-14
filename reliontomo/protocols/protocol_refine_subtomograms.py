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
import glob
import re
from enum import Enum
from typing import Tuple

from pwem.convert.headers import fixVolume
from pwem.objects import SetOfFSCs
from reliontomo.objects import RelionSetOfPseudoSubtomograms
from reliontomo.protocols.protocol_base_refine import ProtRelionRefineBase, IS_RELION_50
from reliontomo import Plugin
from os.path import getmtime
from pyworkflow.protocol import PointerParam, LEVEL_ADVANCED, FloatParam, StringParam, BooleanParam, EnumParam
from reliontomo.constants import ANGULAR_SAMPLING_LIST, SYMMETRY_HELP_MSG, \
    REFINE_FSC_REF_STAR, REFINE_STAR_FSC_TABLE, REFINE_STAR_FSC_COLUMNS
from reliontomo.utils import getProgram
from tomo.objects import AverageSubTomogram


class outputObjects(Enum):
    relionParticles = RelionSetOfPseudoSubtomograms
    average = AverageSubTomogram
    outputFSC = SetOfFSCs


class ProtRelionRefineSubtomograms(ProtRelionRefineBase):
    """
    Performs automated high-resolution subtomogram refinement using
    gold-standard refinement strategies within Relion tomography workflows.
    The protocol iteratively refines particle orientations, positions, and
    three-dimensional reconstructions while minimizing overfitting through
    independent half-set refinement and Fourier Shell Correlation analysis.

    AI Generated:

    3D Auto-Refine for Subtomograms (ProtRelionRefineSubtomograms) —
    User Manual

        Overview

        This protocol performs automated three-dimensional refinement of
        subtomogram datasets using Relion-based gold-standard refinement
        procedures. Its primary purpose is to improve the resolution and
        structural consistency of subtomogram averages by iteratively refining
        particle orientations, translations, and reconstruction parameters in
        a statistically robust manner.

        In cryo-electron tomography workflows, this refinement stage is one
        of the most biologically important steps because it determines the
        final structural quality of the reconstructed macromolecular complex.
        The protocol is designed to operate with minimal user intervention
        while automatically monitoring convergence and avoiding overfitting.

        The refinement strategy relies on independent half-dataset
        reconstructions and Fourier Shell Correlation estimation. This
        gold-standard methodology helps ensure that the reported resolution
        reflects reproducible structural information rather than noise
        amplification. For biological users, this provides greater confidence
        that visible structural features correspond to genuine molecular
        detail.

        General Workflow and Inputs

        The protocol requires a set of subtomogram particles together with an
        initial reference map. The reference provides the starting structural
        framework toward which the particles are iteratively aligned and
        refined. In practice, the initial reference often originates from a
        previous subtomogram averaging stage, a low-resolution reconstruction,
        or a biologically related structure.

        The quality of the starting reference strongly influences refinement
        stability. In most biological workflows, it is advisable to use a
        reference that captures the overall molecular shape without containing
        excessive high-resolution detail that could bias the refinement.

        Since subtomogram refinement is computationally demanding, especially
        for pseudo-subtomogram representations, it is often practical to begin
        refinement at lower resolutions or higher binning levels before
        progressively refining at higher sampling rates. This staged approach
        improves computational efficiency while maintaining stable alignment
        behavior.

        Reference Scaling and Low-Pass Filtering

        The protocol assumes that the reference map and experimental particles
        are represented on compatible intensity scales and spatial sampling.
        Proper scaling is biologically important because inconsistent density
        normalization may interfere with statistical comparison between the
        reference and experimental data.

        When uncertainty exists regarding the intensity normalization of the
        reference map, the protocol can initially rely on cross-correlation
        approaches before transitioning into fully probabilistic refinement.
        This strategy improves robustness when references originate from
        external software packages or differently normalized workflows.

        Low-pass filtering of the initial reference is strongly recommended.
        From a biological perspective, aggressive filtering reduces the risk
        that noise or reference bias dominates early refinement iterations.
        A low-resolution starting reference encourages the refinement to
        recover structural information directly from the experimental data.

        Symmetry Considerations

        The protocol supports the application of molecular symmetry during
        refinement. Correct symmetry selection can substantially improve
        reconstruction quality because symmetry-related particle views are
        combined to increase the effective signal-to-noise ratio.

        For highly symmetric biological assemblies such as viral capsids,
        molecular cages, or oligomeric complexes, applying the appropriate
        symmetry may significantly enhance achievable resolution. However,
        imposing incorrect symmetry can distort biologically meaningful
        structural differences and generate misleading reconstructions.

        Asymmetric refinement should therefore be used whenever structural
        heterogeneity, flexibility, or asymmetry is biologically relevant.
        The protocol also supports symmetry relaxation approaches, which are
        particularly useful for pseudo-symmetric systems where related
        subunits adopt distinct conformations.

        Masking and Solvent Treatment

        Masking plays a central biological role in subtomogram refinement
        because it determines which regions contribute most strongly to the
        alignment and reconstruction process. Proper masking improves signal
        isolation and reduces the influence of surrounding solvent noise.

        A primary solvent mask can be provided to focus refinement on the
        biologically relevant molecular region. In many workflows, soft masks
        surrounding the particle improve refinement stability while minimizing
        sharp boundary artifacts.

        In specialized applications such as viral reconstructions, additional
        masking strategies may be used to flatten internal density regions or
        suppress solvent variability. These approaches can improve refinement
        robustness in highly symmetric particles containing large enclosed
        volumes.

        Solvent-corrected Fourier Shell Correlation calculations may also be
        enabled. These corrections often produce more realistic resolution
        estimates when masks occupy relatively small regions of the
        reconstruction box.

        Angular Sampling and Alignment Strategy

        The protocol employs adaptive angular sampling strategies that
        progressively refine particle orientations during optimization. Early
        iterations typically explore broader orientation ranges, while later
        stages focus on increasingly precise local searches.

        This adaptive refinement is biologically important because it balances
        robustness and computational efficiency. Broad searches improve the
        ability to recover correct particle orientations, whereas fine local
        searches maximize achievable structural detail once convergence is
        approached.

        Advanced users may control local search thresholds and angular
        sampling behavior to optimize difficult datasets. Finer angular
        sampling can accelerate convergence toward high-resolution solutions,
        although excessively aggressive refinement may occasionally reduce
        robustness in heterogeneous samples.

        Translation searches are also refined iteratively. Accurate particle
        centering is essential because even small translational inaccuracies
        can blur high-resolution structural information in the final average.

        Computational Considerations

        Subtomogram refinement is computationally intensive due to the large
        volume of three-dimensional data and the iterative nature of the
        optimization. The protocol supports parallel processing and GPU
        acceleration to improve execution efficiency in large datasets.

        Fourier-space padding options influence interpolation quality during
        refinement. Higher padding generally improves numerical accuracy and
        reduces reconstruction artifacts, although it increases memory usage
        and computational cost.

        In routine biological workflows, default computational parameters are
        often sufficient. However, large particles, high-resolution targets,
        or very large datasets may require careful balancing of computational
        efficiency and reconstruction quality.

        Outputs and Biological Interpretation

        The protocol produces refined subtomogram particles, reconstructed
        three-dimensional maps, independently refined half-maps, and Fourier
        Shell Correlation curves. Together, these outputs provide the basis
        for biological interpretation and downstream structural analysis.

        The final reconstruction represents the consensus structural state
        recovered from the particle population. Improved resolution may reveal
        secondary structure elements, ligand binding regions, membrane
        organization, or conformational variability that were not visible in
        earlier processing stages.

        Half-maps are particularly important because they enable independent
        validation of reconstruction quality and support subsequent
        post-processing procedures such as sharpening and local resolution
        estimation.

        Fourier Shell Correlation curves provide quantitative estimates of
        reconstruction resolution and refinement convergence. Biological users
        should interpret these values together with visual inspection of the
        density map rather than relying exclusively on nominal resolution
        numbers.

        Practical Recommendations

        In most subtomogram averaging workflows, refinement should begin with
        conservative parameters and well-filtered references. Attempting to
        refine aggressively from poor initial models often leads to unstable
        alignments or reference bias.

        Careful masking usually provides one of the largest improvements in
        refinement quality. Masks should include the stable molecular core
        while excluding surrounding solvent regions and highly flexible
        densities whenever possible.

        For challenging datasets containing structural heterogeneity,
        flexibility, or low signal-to-noise ratios, iterative refinement with
        progressively improved references often yields better results than a
        single aggressive refinement run.

        Final Perspective

        Automated subtomogram refinement is not merely a computational
        optimization procedure but a central biological interpretation step in
        cryo-electron tomography. The reliability of the final reconstruction
        depends strongly on appropriate reference preparation, masking,
        symmetry selection, and careful validation of refinement behavior.

        When applied thoughtfully, this protocol enables recovery of highly
        detailed structural information directly within native cellular or
        molecular environments, providing biologically meaningful insight into
        macromolecular organization and function.
    """

    _label = '3D auto-refine'
    _possibleOutputs = outputObjects
    FILE_KEYS = ['data', 'optimiser', 'sampling']
    PREFIXES = ['half1_', 'half2_']

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        self._defineInputParams(form)
        self._defineReferenceParams(form)
        self._defineCTFParams(form)
        self._defineOptimisationParamsCommon2All(form)
        self._defineAutoSamplingParams(form)
        self._defineComputeParams(form)
        self._insertGpuParams(form)
        self._defineAdditionalParams(form)
        form.addParallelSection(threads=0, mpi=3)

    def _defineInputParams(self, form):
        self._defineIOParams(form)
        form.addParam('referenceVolume', PointerParam,
                      pointerClass='Volume',
                      allowsNull=False,
                      label='Reference volume',
                      help='Initial reference 3D map, it should have the same dimensions and the same '
                           'pixel size as your input particles.')
        form.addParam('solventMask', PointerParam,
                      pointerClass='VolumeMask',
                      label='Reference mask (optional)',
                      allowsNull=True,
                      help='A volume mask containing a (soft) mask with the same dimensions '
                           'as the reference(s), and values between 0 and 1, with 1 being 100% protein '
                           'and 0 being 100% solvent. The reconstructed reference map will be multiplied '
                           'by this mask. If no mask is given, a soft spherical mask based on the <radius> '
                           'of the mask for the experimental images will be applied.\n\n'
                           'In some cases, for example for non-empty icosahedral viruses, it is also useful '
                           'to use a second mask. Check _Advaced_ parameters to select another volume mask.')
        if not IS_RELION_50:
            form.addParam('solventMask2', PointerParam,
                          pointerClass='VolumeMask',
                          expertLevel=LEVEL_ADVANCED,
                          allowsNull=True,
                          label='Second reference mask (optional)',
                          help='For all white (value 1) pixels in this second mask the '
                               'corresponding pixels in the reconstructed map are set to the average value of '
                               'these pixels. Thereby, for example, the higher density inside the virion may be '
                               'set to a constant. Note that this second mask should have one-values inside the '
                               'virion and zero-values in the capsid and the solvent areas.')

    def _defineReferenceParams(self, form):
        form.addSection(label='Reference')
        form.addParam('isMapAbsoluteGreyScale', BooleanParam,
                      default=True,
                      label='Is initial 3D map on absolute greyscale?',
                      help='Probabilities are calculated based on a Gaussian noise model,'
                           'which contains a squared difference term between the reference and the experimental '
                           'image.\n\n This has a consequence that the reference needs to be on the same absolute '
                           'intensity greyscale as the experimental images. RELION and XMIPP reconstruct maps at '
                           'their absolute intensity greyscale. Other packages may perform internal normalisations of '
                           'the reference density, which will result in incorrect grey-scales. But, if the map was'
                           'reconstructed in RELION or in XMIPP, set this option to Yes, otherwise set it to No.\n\n'
                           'If set to No, RELION will use a (grey-scale invariant) cross-correlation criterion in the '
                           'first iteration, and prior to the second iteration the map will be filtered again using '
                           'the initial low-pass filter. This procedure is relatively quick and typically does not '
                           'negatively affect the outcome of the subsequent map refinement. Therefore, if in doubt it '
                           'is recommended to set this option to No.')
        if IS_RELION_50:
            form.addParam('doResizeRef', BooleanParam,
                          default=True,
                          label='Resize references if needed?',
                          help='If true, and if the input reference map (and mask) do not have the same pixel size '
                               'and/or box size, then they will be re-scaled and re-boxed accordingly. If this option '
                               'is set to false, then the program will die with an error if the reference does not '
                               'have the correct pixel and/or box size.')
        form.addParam('initialLowPassFilterA', FloatParam,
                      default=30,
                      label='Initial low-pass filter (A)',
                      help='It is recommended to strongly low-pass filter your initial reference map. '
                           'If it has not yet been low-pass filtered, it may be done internally using this option. '
                           'If set to 0, no low-pass filter will be applied to the initial reference(s).')

        help3drefine = 'If the molecule is asymmetric, set Symmetry group to C1. Note their are multiple possibilities' \
                       ' for icosahedral symmetry: \n ' \
                       ' _* I1_: No-Crowther 222 (standard in Heymann, Chagoyen & Belnap, JSB, 151 (2005) (196-207)\n' \
                       ' _* I2_: Crowther 222\n ' \
                       ' _* I3_: 52-setting (as used in SPIDER?) \n' \
                       ' _* I4_: A different 52 setting \n' \
                       'RELION uses XMIPPs libraries for symmetry operations. Therefore, look at the XMIPP:\n'+ SYMMETRY_HELP_MSG
        super()._insertSymmetryParam(form, help3drefine)

    def _defineOptimisationParamsCommon2All(self, form):
        self._insertOptimisationSection(form)
        self._insertMaskDiameterParam(form)
        self._insertZeroMaskParam(form)
        form.addParam('solventCorrectFSC', BooleanParam,
                      default=False,
                      condition='solventMask',
                      label='Use solvent-flattened FSCs?',
                      help="If set to Yes, then instead of using unmasked maps to calculate the gold-standard FSCs "
                           "during refinement, masked half-maps are used and a post-processing-like correction of "
                           "the FSC curves (with phase-randomisation) is performed every iteration.\n\n"
                           "This only works when a reference mask is provided. This may yield "
                           "higher-resolution maps, especially when the mask contains only a relatively small "
                           "volume inside the box.")
        if IS_RELION_50:
            self._insertBlushRegParam(form)

    def _defineAutoSamplingParams(self, form):
        form.addSection(label='Auto-sampling')
        super()._insertAngularCommonParams(form)
        form.addParam('localSearchAutoSampling', EnumParam,
                      default=4,
                      choices=ANGULAR_SAMPLING_LIST,
                      label='Local searches from auto-sampling',
                      help="Minimum healpix order (before oversampling) from which autosampling procedure will "
                           "use local searches.\n\n"
                           "In the automated procedure to increase the angular samplings, local angular "
                           "searches of -6/+6 times the sampling rate will be used from this angular sampling rate "
                           "onwards. For most lower-symmetric particles a value of 1.8 degrees will be sufficient. "
                           "Perhaps icosahedral symmetries may benefit from a smaller value such as 0.9 degrees.")
        self._insertRelaxSymmetry(form)
        form.addParam('useFinerAngularSampling', BooleanParam,
                      default=False,
                      label='Use finer angular sampling faster?',
                      help="If set to Yes, then let auto-refinement proceed faster with finer angular samplings. "
                           "Two additional conditions will be considered:\n\n "
                           "\t-Angular sampling will go down despite changes still happening in the angles. (this is "
                           "the Relion flag  --auto_ignore_angles )\n"
                           "\t-Angular sampling will go down if the current resolution already requires that sampling "
                           "(this is the Relion flag --auto_resol_angles)\n"
                           "\t at the edge of the particle.\n\nThis option will make the computation faster, but "
                           "hasn't been tested for many cases for potential loss in reconstruction quality upon "
                           "convergence.")
        if IS_RELION_50:
            self._insertPriorWidthParam(form)

    def _defineComputeParams(self, form, isOnlyClassif=False):
        super()._defineComputeParams(form, isOnlyClassif=isOnlyClassif)
        form.addParam('skipPadding', BooleanParam,
                      default=False,
                      label='Skip padding?',
                      help="If set to Yes, the calculations will not use padding in Fourier space for better "
                           "interpolation in the references. Otherwise, references are padded 2x before Fourier "
                           "transforms are calculated. Skipping padding (i.e. use --pad 1) gives nearly as good "
                           "results as using --pad 2, but some artifacts may appear in the corners from signal "
                           "that is folded back.")

    @staticmethod
    def _insertRelaxSymmetry(form, condition=True):
        form.addParam('relaxSym', StringParam,
                      allowsNull=True,
                      condition=condition,
                      label='Symmetry to be relaxed',
                      help="With this option, poses related to the standard local angular search range by the given "
                           "point group will also be explored. For example, if you have a pseudo-symmetric dimer A-A', "
                           "refinement or classification in C1 with symmetry relaxation by C2 might be able to improve "
                           "distinction between A and A'. Note that the reference must be more-or-less aligned to the "
                           "convention of (pseudo-)symmetry operators. For details, see Ilca et al 2019 and Abrishami "
                           "et al 2020 cited in the About dialog.\n\n")

    # -------------------------- INSERT steps functions -----------------------
    def _insertAllSteps(self):
        self._initialize()
        self._insertFunctionStep(self.convertInputStep, needsGPU=False)
        self._insertFunctionStep(self.autoRefineStep, needsGPU=True)
        self._insertFunctionStep(self.createOutputStep, needsGPU=False)

    # -------------------------- STEPS functions ------------------------------
    def _initialize(self):
        """ This function is meant to be called after the
        working dir for the protocol have been set. (maybe after recovery from mapper)
        """
        self._createFilenameTemplates()
        self._createIterTemplates()

    def autoRefineStep(self):
        nMpi = self.numberOfMpi.get()
        Plugin.runRelionTomo(self, getProgram('relion_refine', nMpi),
                             self._genAutoRefineCommand(),
                             numberOfMpi=nMpi)

    def createOutputStep(self):
        inParticles = self.getInputParticles()

        # Output Relion particles
        relionParticles = self.genRelionParticles(optimisationFileName='_optimisation_set.star',
                                                  particles='_data.star')

        # Output volume
        vol = AverageSubTomogram()
        volName = self._getRefineResultFn()
        fixVolume(volName)  # Fix header for xmipp to consider it a volume instead of a stack
        vol.setFileName(volName)
        vol.setSamplingRate(relionParticles.getCurrentSamplingRate())
        half1, half2 = self._getRefineResultHalves()
        vol.setHalfMaps([half1, half2])

        # Output FSC
        fn = self._getExtraPath(REFINE_FSC_REF_STAR)

        setOfFSC = self.genFSCs(fn, REFINE_STAR_FSC_TABLE,
                                REFINE_STAR_FSC_COLUMNS)

        outputDict = {outputObjects.relionParticles.name: relionParticles,
                      outputObjects.average.name: vol,
                      outputObjects.outputFSC.name: setOfFSC}
        self._defineOutputs(**outputDict)
        self._defineSourceRelation(inParticles, relionParticles)
        self._defineSourceRelation(inParticles, vol)
        self._defineSourceRelation(inParticles, setOfFSC)

    # -------------------------- INFO functions -------------------------------

    # --------------------------- UTILS functions -----------------------------
    def _getRefineResultFn(self) -> str:
        return self._getExtraPath('_class001.mrc')

    def _getRefineResultHalves(self) -> Tuple[str, str]:
        pattern = '*it*half%s_class*.mrc'
        half1 = self._getLastFileName(self._getExtraPath(pattern % 1))
        half2 = self._getLastFileName(self._getExtraPath(pattern % 2))
        return half1, half2

    def _genAutoRefineCommand(self):
        cmd = self._genBaseCommand()
        cmd += ' --auto_refine --split_random_halves --low_resol_join_halves 40 --norm --scale --flatten_solvent '
        # I/O args

        ref = self.referenceVolume.get()
        refFile = ref.getFileName()
        # If halves exists, use a half as explained here:
        # https://relion.readthedocs.io/en/release-4.0/STA_tutorial/Refine3D.html#high-resolution-3d-refinement
        if ref.hasHalfMaps():
            refFile = ref._halfMapFilenames[0]

        cmd += '--ref %s ' % refFile
        if self.solventMask.get():
            cmd += '--solvent_mask %s ' % self.solventMask.get().getFileName()
        # Reference args
        if not self.isMapAbsoluteGreyScale.get():
            cmd += '--firstiter_cc '
        if self.initialLowPassFilterA.get():
            cmd += '--ini_high %.2f ' % self.initialLowPassFilterA.get()
        cmd += '--sym %s ' % self.symmetry.get()
        # Optimisation args
        if self.solventCorrectFSC.get():
            cmd += '--solvent_correct_fsc '
        if self.zeroMask.get():
            cmd += '--zero_mask '
        # Angular sampling args
        cmd += '--healpix_order %i ' % self.angularSamplingDeg.get()  # - self.oversampling.get())
        cmd += '--offset_range %i ' % self.offsetSearchRangePix.get()
        cmd += '--offset_step %i ' % (self.offsetSearchStepPix.get() * 2 ** self.oversampling.get())
        cmd += '--auto_local_healpix_order %i ' % self.localSearchAutoSampling.get()  # - self.oversampling.get())
        if self.relaxSym.get():
            cmd += '--relax_sym %s ' % self.relaxSym.get()
        if self.useFinerAngularSampling.get():
            cmd += '--auto_ignore_angles --auto_resol_angles '
        # Compute args
        cmd += '--pad %i ' % (1 if self.skipPadding.get() else 2)
        # Plugin version-dependent parameters
        if IS_RELION_50:
            # Auto-sampling
            cmd += '--sigma_tilt %i ' % self.priorWidthTiltAngle.get()
            # Reference
            if self.doResizeRef.get():
                cmd += '--trust_ref_size '
            # Optimisation
            if self.doBlushReg.get():
                cmd += '--blush '
        else:
            # Input
            if self.solventMask2.get():
                cmd += '--solvent_mask2 %s ' % self.solventMask2.get().getFileName()
                # Compute
                if not self.skipGridding.get():
                    cmd += '--dont_skip_gridding '

        return cmd

    @staticmethod
    def _getLastFileName(pattern):
        files = glob.glob(pattern)
        files.sort(key=getmtime)
        return files[-1]

    def _createFilenameTemplates(self):
        """ Centralize how files are called for iterations and references. """
        self.extraIter = self._getExtraPath('_it%(iter)03d_')
        myDict = {
            'input_star': self._getPath('input_particles.star'),
            'input_mrcs': self._getPath('input_particles.mrcs'),
            'data_scipion': self.extraIter + 'data_scipion.sqlite',
            'projections': self.extraIter + '%(half)sclass%(ref3d)03d_projections.sqlite',
            'classes_scipion': self.extraIter + 'classes_scipion.sqlite',
            'data': self.extraIter + 'data.star',
            'model': self.extraIter + 'model.star',
            'optimiser': self.extraIter + 'optimiser.star',
            'angularDist_xmipp': self.extraIter + 'angularDist_xmipp.xmd',
            'all_avgPmax': self._getPath('iterations_avgPmax.star'),
            'all_changes': self._getPath('iterations_changes.star'),
            'selected_volumes': self._getPath('selected_volumes_xmipp.xmd'),
            'dataFinal': self._getExtraPath("_data.star"),
            'modelFinal': self._getExtraPath("_model.star"),
            'finalvolume': self._getExtraPath("_class%(ref3d)03d.mrc"),
            'final_half1_volume': self._getExtraPath("_half1_class%(ref3d)03d_unfil.mrc"),
            'final_half2_volume': self._getExtraPath("_half2_class%(ref3d)03d_unfil.mrc"),
            'finalSGDvolume': self._getExtraPath("r_it%(iter)03d_class%(ref3d)03d.mrc"),
            'preprocess_particles': self._getPath("preprocess_particles.mrcs"),
            'preprocess_particles_star': self._getPath("preprocess_particles.star"),
            'preprocess_particles_prefix': "preprocess_particles"
        }
        # add to keys, data.star, optimiser.star and sampling.star
        for key in self.FILE_KEYS:
            myDict[key] = self.extraIter + '%s.star' % key
            key_xmipp = key + '_xmipp'
            myDict[key_xmipp] = self.extraIter + '%s.xmd' % key
        # add other keys that depends on prefixes
        for p in self.PREFIXES:
            myDict['%smodel' % p] = self.extraIter + '%smodel.star' % p
            myDict['%svolume' % p] = self.extraIter + p + 'class%(ref3d)03d.mrc'

        self._updateFilenamesDict(myDict)

    def _createIterTemplates(self):
        """ Setup the regex on how to find iterations. """
        self._iterTemplate = self._getFileName('data', iter=0).replace('000', '???')
        # Iterations will be identify by _itXXX_ where XXX is the iteration number
        # and is restricted to only 3 digits.
        self._iterRegex = re.compile('_it(\d{3})_')

    # -------------------------- INFO functions -------------------------------
