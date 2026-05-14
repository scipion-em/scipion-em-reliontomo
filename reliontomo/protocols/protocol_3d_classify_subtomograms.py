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
from os import remove, listdir
from os.path import abspath, exists, isfile, join

from emtable import Table

from pyworkflow.protocol import Form
from pyworkflow.object import Float
from pyworkflow.utils import moveFile, createLink
from reliontomo.constants import OUT_PARTICLES_STAR, TOMO_PARTICLE_ID, OPTICS_TABLE, PARTICLES_TABLE
from reliontomo.objects import RelionSetOfPseudoSubtomograms
from reliontomo.protocols import ProtRelionRefineSubtomograms
from reliontomo.protocols.protocol_base_refine import ProtRelionRefineBase
from reliontomo import Plugin
from pyworkflow.protocol import FloatParam, BooleanParam, GE, LE, IntParam, StringParam
from reliontomo.protocols.protocol_base_relion import IS_RELION_50
from reliontomo.utils import getProgram
from tomo.objects import SetOfClassesSubTomograms, SetOfSubTomograms


class outputObjects(Enum):
    relionParticles = RelionSetOfPseudoSubtomograms
    classes = SetOfClassesSubTomograms


class ProtRelion3DClassifySubtomograms(ProtRelionRefineSubtomograms):
    """
    Performs 3D classification of subtomograms in cryo-electron tomography workflows. The protocol separates
    heterogeneous particle populations into distinct structural classes while optionally refining particle
    orientations and translations during the classification process.

    AI Generated:

    3D Classification of Subtomograms (ProtRelion3DClassifySubtomograms) - User Manual
        Overview

        The 3D Classification of Subtomograms protocol is designed to identify structural variability within
        cryo-electron tomography datasets by grouping subtomograms into multiple biologically meaningful
        classes. This process is particularly important when the sample contains conformational flexibility,
        compositional heterogeneity, partially assembled complexes, or particles in different functional states.

        In practical cryo-ET studies, classification is often the key step that separates a mixed population
        into cleaner and more interpretable subsets. Instead of producing a single consensus reconstruction,
        the protocol attempts to reveal the underlying structural diversity present in the dataset. This is
        especially valuable for studying dynamic molecular machines, membrane-associated complexes, or crowded
        cellular environments where particles rarely exist in a single homogeneous state.

        Inputs and General Workflow

        The protocol requires a set of subtomograms together with an initial reference volume. The reference
        provides the starting structural information used during classification and alignment. In many cases,
        this reference originates from a previous subtomogram refinement, an averaged reconstruction, or an
        external structure filtered to an appropriate resolution.

        During processing, the protocol iteratively compares each particle against several evolving class
        references. Over successive iterations, particles are assigned to the classes that best represent
        their structural features. As the refinement progresses, the class averages gradually diverge and
        become more specialized, allowing biologically distinct conformations or assemblies to emerge.

        The number of requested classes is biologically important. Using too few classes may merge distinct
        structural states into a blurred average, while using too many classes may fragment the dataset into
        unstable or poorly populated groups. In exploratory studies, it is common to test several class
        numbers before selecting the most biologically meaningful solution.

        Classification Versus Alignment

        The protocol can operate either with simultaneous alignment and classification or with classification
        alone. Simultaneous alignment is generally recommended when particle orientations are still uncertain.
        In this mode, the protocol refines orientations and translations while separating particles into
        classes, improving the consistency of each resulting reconstruction.

        Classification without alignment is useful when reliable particle orientations already exist from a
        previous refinement. This strategy is particularly effective when using highly focused masks to study
        local variability within a stable complex. For example, a user may classify only the flexible domain
        of a ribosome, membrane receptor, or viral spike while preserving previously determined global
        orientations.

        Angular Sampling and Local Searches

        Angular sampling parameters control how broadly the protocol explores particle orientations during
        refinement. Wider searches provide greater robustness when orientations are uncertain, although they
        increase computational cost. Narrower searches are more efficient when particles are already roughly
        aligned.

        Local angular searches are especially useful during late refinement stages or when studying subtle
        conformational differences. Instead of exploring the full orientation space, the protocol restricts
        the search around previously determined orientations. This approach improves stability and reduces
        computational demands while preserving biologically relevant refinements.

        The protocol also supports adaptive sampling strategies that begin with coarser searches and become
        progressively more precise as the refinement converges. This behavior is particularly advantageous
        for large datasets where computational efficiency is important.

        Reference Masks and Solvent Treatment

        Masking plays a central role in obtaining biologically meaningful classes. A solvent mask allows the
        protocol to focus classification on structurally relevant regions while suppressing noise from solvent
        areas or flexible background density. Soft masks are generally preferred because they reduce edge
        artifacts and improve refinement stability.

        Focused masking becomes especially important when studying local heterogeneity. In many biological
        systems, only a small region of the structure changes between states, while the remainder stays
        relatively constant. By concentrating classification on the variable region, the protocol becomes
        more sensitive to subtle conformational differences.

        Solvent flattening and particle masking options further improve the stability of the refinement by
        reducing the influence of noisy regions outside the particle boundary. Depending on the dataset,
        either zero-filled or noise-filled solvent regions may provide better performance.

        Regularisation and Overfitting Control

        Classification protocols are inherently more susceptible to overfitting than gold-standard refinement
        procedures. For this reason, regularisation parameters and resolution limits are biologically
        important safeguards. Increasing regularisation generally favors smoother and more stable classes,
        whereas weaker regularisation may capture finer details at the risk of amplifying noise.

        The protocol allows limiting the resolution used during alignment calculations. This option is often
        beneficial for difficult datasets with small particles, low contrast, or high flexibility. Restricting
        the alignment to moderate resolutions can stabilize classification and prevent the algorithm from
        fitting noise instead of genuine structural features.

        Neural-network-based regularisation methods may also be available in modern processing environments.
        These approaches attempt to denoise intermediate reconstructions during refinement and can improve
        class interpretability for challenging datasets.

        Computational Considerations

        Subtomogram classification is computationally demanding because each iteration involves repeated
        comparisons between particles and multiple class references. GPU acceleration can significantly reduce
        runtime and is strongly recommended for large datasets.

        The protocol supports several performance optimization strategies, including pooled particle handling,
        adaptive oversampling, parallel execution, and accelerated subset-based refinement. Fast subset modes
        are particularly useful during early exploratory analysis because they allow rapid estimation of the
        main structural variability before committing to a full refinement.

        Efficient storage systems and adequate memory resources are also important because tomography datasets
        can become extremely large, especially when using many particles or large box sizes.

        Outputs and Their Interpretation

        The protocol produces a refined particle set together with a collection of 3D classes representing
        the structural populations identified during processing. Each class includes a representative volume
        and associated metadata describing particle assignments and refinement statistics.

        Biologically, the resulting classes should be interpreted carefully. Well-resolved and highly populated
        classes often correspond to stable conformational or compositional states. Smaller or noisier classes
        may represent rare intermediates, damaged particles, misaligned subsets, or residual heterogeneity.

        Class populations should not automatically be interpreted as exact biological abundances because
        classification quality can depend on masking, alignment stability, particle quality, and refinement
        parameters. Visual inspection and downstream biological validation remain essential.

        Practical Recommendations

        In most biological workflows, it is advisable to begin with a moderate number of classes and relatively
        conservative alignment settings. Excessively aggressive classification may fragment the dataset and
        complicate interpretation. Once the major conformational states become visible, additional focused
        classifications can be performed to study specific regions in greater detail.

        For highly heterogeneous samples, applying an appropriate solvent mask often provides the largest
        improvement in class quality. Local angular searches are particularly effective when the dataset has
        already undergone a reliable global refinement.

        When classification results appear unstable or noisy, reducing alignment resolution, increasing
        regularisation, or simplifying the number of requested classes often improves robustness. Conversely,
        well-behaved datasets with strong signal may benefit from finer angular sampling and more detailed
        local searches.

        Final Perspective

        Subtomogram classification is one of the most biologically informative stages of cryo-electron
        tomography analysis because it transforms heterogeneous particle populations into interpretable
        structural states. Careful selection of masks, alignment strategies, regularisation parameters,
        and class numbers is essential for revealing meaningful biological variability while minimizing
        overfitting and noise amplification.
    """

    _label = '3D classification'
    modelTable = Table()
    classesTable = Table()
    opticsTable = Table()
    particlesTable = Table()
    _possibleOutputs = outputObjects

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        self._defineInputParams(form)
        self._insertBinThreadsParam(form)
        self._defineReferenceParams(form)
        self._defineCTFParams(form)
        self._defineOptimisationParams(form)
        self._defineSamplingParams(form)
        self._defineComputeParams(form)
        self._insertGpuParams(form)
        self._defineAdditionalParams(form)
        form.addParallelSection(threads=0, mpi=3)

    def _defineOptimisationParams(self, form):
        self._insertOptimisationSection(form)
        self._insertNumOfClassesParam(form)
        self._insertRegularisationParam(form, isCl3d=True)
        form.addParam('nIterations', IntParam,
                      label='Number of iterations',
                      default=25,
                      help='Number of iterations to be performed. Note that the current implementation of 2D class '
                           'averaging and 3D classification does NOT comprise a convergence criterium. Therefore, '
                           'the calculations will need to be stopped by the user if further iterations do not yield '
                           'improvements in resolution or classes.')
        form.addParam('useFastSubsets', BooleanParam,
                      label='Use fast subsets (for large data sets)?',
                      default=False,
                      help='If set to Yes, the first 5 iterations will be done with random subsets of only K*1500 '
                           'particles (K being the number of classes); the next 5 with K*4500 particles, the next '
                           '5 with 30% of the data set; and the final ones with all data. This was inspired by a '
                           'cisTEM implementation by Niko Grigorieff et al.')
        self._insertMaskDiameterParam(form)
        self._insertZeroMaskParam(form)
        form.addParam('limitResolutionEStepTo', FloatParam,
                      label='Limit resolution E-step to (Å)',
                      default=-1,
                      help='If set to a positive number, then the expectation step (i.e. the alignment) will be done '
                           'only including the Fourier components up to this resolution (in Angstroms).\n\nThis is '
                           'useful to prevent overfitting, as the classification runs in RELION are not guaranteed to '
                           'be 100% overfitting-free (unlike the 3D auto-refine with its gold-standard FSC). '
                           'In particular for very difficult data sets, e.g. of very small or featureless particles, '
                           'this has been shown to give much better class averages. \n\nIn such cases, values in the '
                           'range of 7-12 Angstroms have proven useful.')
        if IS_RELION_50:
            self._insertBlushRegParam(form)

    def _defineSamplingParams(self, form: Form):
        form.addSection(label='Sampling')
        form.addParam('doImageAlignment', BooleanParam,
                      label='Perform image alignment?',
                      default=False,
                      help='If set to No, then rather than performing both alignment and classification, only '
                           'classification will be performed. This allows the use of very focused masks. It requires '
                           'that the optimal orientations of all particles are already stored in the input STAR file.')

        self._insertAngularCommonParams(form, condition='doImageAlignment')

        form.addParam('doLocalAngleSearch', BooleanParam,
                      label='Perform local angular searches?',
                      default=False,
                      condition='doImageAlignment',
                      help="If set to Yes, then rather than performing exhaustive angular searches, local searches "
                           "within the range given below will be performed.\n\nA prior Gaussian distribution centered "
                           "at the optimal orientation in the previous iteration and with a stddev of 1/3 of the range "
                           "given below will be enforced.")
        form.addParam('localAngularSearchRange', FloatParam,
                      label='Local angular search range (deg.)',
                      condition='doImageAlignment and doLocalAngleSearch',
                      default=5,
                      validators=[GE(0), LE(15)],
                      help="Local angular searches will be performed within +/- the given amount (in degrees) from "
                           "the optimal orientation in the previous iteration.\n\nA Gaussian prior (also see previous "
                           "option) will be applied, so that orientations closer to the optimal orientation in the "
                           "previous iteration will get higher weights than those further away.")

        self._insertRelaxSymmetry(form, condition='doImageAlignment and doLocalAngleSearch')
        form.addParam('allowCoarser', BooleanParam,
                      label='Allow coarser sampling?',
                      default=False,
                      condition='doImageAlignment',
                      help="If set to Yes, the program will use coarser angular and translational samplings if the "
                           "estimated accuracies of the assignments are still low in the earlier iterations. This may "
                           "speed up the calculations.")
        if IS_RELION_50:
            self._insertPriorWidthParam(form)

    # -------------------------- INSERT steps functions -----------------------
    def _insertAllSteps(self):
        self._insertFunctionStep(self.convertInputStep, needsGPU=False)
        self._insertFunctionStep(self._classify3d, needsGPU=True)
        self._insertFunctionStep(self.createOutputStep, needsGPU=False)

    # -------------------------- STEPS functions ------------------------------
    def _classify3d(self):
        nMpi = self.numberOfMpi.get()
        Plugin.runRelionTomo(self, getProgram('relion_refine', nMpi), self._genCl3dCommand(), numberOfMpi=nMpi)

    def createOutputStep(self):
        inParticles = self.inReParticles.get()

        # Rename the particles file generated (_data.star) to follow the name convention
        createLink(self._getIterGenFileName('data', self.nIterations.get()), self._getExtraPath(OUT_PARTICLES_STAR))

        # Output RelionParticles
        relionParticles = self.genRelionParticles()

        # Output classes
        classes3D = self._createSetOfClassesSubTomograms(relionParticles)
        classes3D.setImages(relionParticles)
        self._fillClassesFromIter(classes3D, self.nIterations.get())

        # Register the outputs
        outputDict = {outputObjects.relionParticles.name: relionParticles,
                      outputObjects.classes.name: classes3D}
        self._defineOutputs(**outputDict)
        self._defineSourceRelation(inParticles, relionParticles)
        self._defineSourceRelation(inParticles, classes3D)

        # Remove undesired files if requested
        if self.keepOnlyLastIterFiles.get():
            self._cleanUndesiredFiles()

    # -------------------------- INFO functions -------------------------------

    # --------------------------- UTILS functions -----------------------------
    def _genCl3dCommand(self):
        cmd = '--norm --scale --flatten_solvent '

        # I/O args
        cmd += self._genIOBaseCmd()
        cmd += '--ref %s ' % self.referenceVolume.get().getFileName()
        if self.solventMask.get():
            cmd += '--solvent_mask %s ' % self.solventMask.get().getFileName()
        if not IS_RELION_50 and self.solventMask2.get():
            cmd += '--solvent_mask2 %s ' % self.solventMask2.get().getFileName()

        # Reference args
        if not self.isMapAbsoluteGreyScale.get():
            cmd += '--firstiter_cc '
        if self.initialLowPassFilterA.get():
            cmd += '--ini_high %.2f ' % self.initialLowPassFilterA.get()
        cmd += '--sym %s ' % self.symmetry.get()
        if IS_RELION_50 and self.doResizeRef.get():
            cmd += '--trust_ref_size '

        # CTF args
        cmd += self._genCTFBaseCmd()

        # Optimisation args
        cmd += self._genOptimisationBaseCmd()
        if self.zeroMask.get():
            cmd += '--zero_mask '
        cmd += '--K %i ' % self.numberOfClasses.get()
        cmd += '--tau2_fudge %d ' % self.regularisation.get()
        cmd += '--iter %i ' % self.nIterations.get()
        if self.useFastSubsets.get():
            cmd += '--fast_subsets '
        if self.limitResolutionEStepTo.get() > 0:
            cmd += '--strict_highres_exp %d ' % self.limitResolutionEStepTo.get()
        if IS_RELION_50 and self.doBlushReg.get():
            cmd += '--blush '

        # Sampling args
        if self.doImageAlignment.get():
            cmd += '--healpix_order %i ' % self.angularSamplingDeg.get()
            cmd += '--offset_range %i ' % self.offsetSearchRangePix.get()
            cmd += '--offset_step %d ' % (self.offsetSearchStepPix.get() * 2 ** self.oversampling.get())
            if self.doLocalAngleSearch.get():
                cmd += '--sigma_ang %d ' % (self.localAngularSearchRange.get() / 3)

                if self.relaxSym.get():
                    cmd += '--relax_sym %s ' % self.relaxSym.get()
            if self.allowCoarser.get():
                cmd += '--allow_coarser_sampling '
            if IS_RELION_50:
                cmd += '--sigma_tilt %i ' % self.priorWidthTiltAngle.get()
        else:
            cmd += '--skip_align '

        # Compute args
        cmd += self._genComputeBaseCmd(onlyCl3d=not self.doImageAlignment.get())
        cmd += '--pad %i ' % (1 if self.skipPadding.get() else 2)

        # Additional args
        cmd += self._genAdditionalBaseCmd()

        return cmd

    def _createSetOfClassesSubTomograms(self, subTomograms: SetOfSubTomograms, suffix=''):
        classes = SetOfClassesSubTomograms.create(self._getPath(),
                                                  template='subtomogramClasses%s.sqlite',
                                                  suffix=suffix)
        classes.setImages(subTomograms)
        classes.setCoordinates3D(subTomograms.getCoordinates3D(asPointer=True))
        return classes

    def _fillClassesFromIter(self, clsSet, iteration):
        """ Create the SetOfClasses3D from a given iteration. """
        self._loadClassifyInfo(iteration)
        clsSet.classifyItems(updateItemCallback=self._updateParticle,
                             updateClassCallback=self._updateClass,
                             itemDataIterator=self.particlesTable.__iter__())

    def _loadClassifyInfo(self, iteration):
        """ Read some information about the produced Relion 3D classes
        from the *model.star file.
        """
        self._classesInfo = {}  # store classes info, indexed by class id
        modelStar = self._getIterGenFileName('model', iteration)
        with open(modelStar) as fid:
            self.modelTable.readStar(fid, tableName='model_general')
            self.classesTable.readStar(fid, tableName='model_classes')
        dataStar = self._getIterGenFileName('data', iteration)
        with open(dataStar) as fid:
            self.opticsTable.readStar(fid, tableName=OPTICS_TABLE)
            self.particlesTable.readStar(fid, tableName=PARTICLES_TABLE)
            if not IS_RELION_50:
                self.particlesTable.sort(TOMO_PARTICLE_ID)

        # Model table has only one row, while classes table has the same number of rows as classes found
        self.nClasses = int(self.modelTable._rows[0].rlnNrClasses)  # self.numberOfClasses.get()
        # Adapt data to the format expected by classifyItems and its callbacks
        for i, row in zip(range(self.nClasses), self.classesTable.__iter__()):
            self._classesInfo[i + 1] = (i + 1, row)

    @staticmethod
    def _updateParticle(item, row):
        item.setClassId(row.rlnClassNumber)  # rlnGroupNumber))
        item._rlnLogLikeliContribution = Float(row.rlnLogLikeliContribution)
        item._rlnMaxValueProbDistribution = Float(row.rlnMaxValueProbDistribution)

    def _updateClass(self, item):
        classId = item.getObjId()
        if classId in self._classesInfo:
            _, row = self._classesInfo[classId]
            fn = row.rlnReferenceImage
            item.setAlignment3D()

            # Representative stuff
            representative = item.getRepresentative()
            representative.setLocation(fn)
            representative.setSamplingRate(self.getInputParticles().getCurrentSamplingRate())
            # relion mrc are technically stacks. Fix this
            representative.fixMRCVolume()

            item._rlnclassDistribution = Float(row.rlnClassDistribution)
            item._rlnAccuracyRotations = Float(row.rlnAccuracyRotations)
            item._rlnAccuracyTranslations = Float(row.rlnAccuracyTranslationsAngst)

    def _cleanUndesiredFiles(self):
        """Remove all files generated by relion_classify 3d excepting the ones which
        correspond to the last iteration. Example for iteration 25:
        it025_class002.mrc
        it025_class001.mrc
        it025_model.star
        it025_sampling.star
        it025_optimiser.star
        it025_data.star
        """
        itPref = 'it'
        clPref = 'class'
        starExt = '.star'
        mrcExt = '.mrc'
        # Classify calculations related files
        calcFiles = ['data', 'model', 'optimiser', 'sampling']
        for i in range(self.nIterations.get()):
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

    def _manageGeneratedFiles(self):
        """There's some kind of bug in relion4 which makes it generate the file in the protocol base directory
        instead of the extra directory. It uses extra as a prefix of each generated file instead. Hence, until
        it's solved, the files will be moved to the extra directory and the prefix extra_ will be removed"""
        prefix = 'extra_'
        genFiles = [f for f in listdir(self._getPath()) if isfile(join(self._getPath(), f))]
        for f in genFiles:
            if f.startswith(prefix):
                moveFile(self._getPath(f), self._getExtraPath(f.replace(prefix, '')))

    def _getIterGenFileName(self, fileType, iteration):
        # Pattern to reproduce is it[3 digit number]_[file type].star
        return self._getExtraPath('_it%03d_%s.star' % (iteration, fileType))
