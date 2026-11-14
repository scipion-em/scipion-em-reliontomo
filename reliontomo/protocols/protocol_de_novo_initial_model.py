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
from os.path import exists
from typing import Union
from pwem.convert.headers import fixVolume
from pwem.objects import Volume, SetOfVolumes
from pyworkflow.object import Pointer
from reliontomo.constants import SYMMETRY_HELP_MSG
from reliontomo.protocols.protocol_base_refine import ProtRelionRefineBase
from reliontomo.protocols.protocol_base_relion import IS_RELION_50
from reliontomo import Plugin
from pyworkflow.protocol import LEVEL_ADVANCED
from reliontomo.utils import getProgram


class outputObjects(Enum):
    average = Volume
    averages = SetOfVolumes


class ProtRelionDeNovoInitialModel(ProtRelionRefineBase):
    """
    Generates a de novo 3D initial model from pseudo-subtomograms in
    cryo-electron tomography workflows. The protocol is designed to
    create an initial low-resolution structural reference without the
    need for a previously known map, enabling downstream refinement,
    classification, and structural interpretation.

    AI Generated:

    De Novo Initial Model (ProtRelionDeNovoInitialModel) — User Manual
        Overview

        The De Novo Initial Model protocol generates one or more initial
        three-dimensional maps directly from pseudo-subtomograms using
        Relion tomography workflows. Its main purpose is to provide an
        unbiased starting structure for subtomogram averaging and
        refinement when no reliable reference model is available.

        In cryo-electron tomography, obtaining an initial reference is
        often one of the most challenging stages of analysis. Biological
        structures may be poorly characterized, structurally flexible,
        or entirely unknown. This protocol addresses that problem by
        reconstructing a low-resolution model directly from the particle
        data itself, allowing users to begin iterative refinement and
        classification procedures without introducing strong external
        bias.

        The protocol is particularly valuable during exploratory studies
        of macromolecular complexes, membrane assemblies, viral particles,
        or heterogeneous cellular structures where structural information
        is incomplete or unavailable.

        Inputs and General Workflow

        The protocol requires a set of pseudo-subtomograms generated from
        tomographic data. These particles should already represent the
        biological structure of interest and should ideally contain a
        broad distribution of orientations. Angular diversity is one of
        the key requirements for successful de novo reconstruction because
        insufficient orientation coverage can lead to distorted or unstable
        maps.

        The workflow begins by estimating an initial low-resolution volume
        directly from the particle set. This model is iteratively improved
        while simultaneously refining particle orientations and translations.
        The resulting structure can then serve as the starting point for
        downstream subtomogram averaging, classification, or high-resolution
        refinement.

        From a biological perspective, the protocol is intended to recover
        the dominant structural organization present in the dataset rather
        than detailed atomic features. The resulting maps should therefore
        be interpreted as low-resolution consensus representations suitable
        for further refinement.

        Symmetry Considerations

        Symmetry handling is one of the most biologically important aspects
        of de novo model generation. Many macromolecular assemblies exhibit
        rotational or point-group symmetry, and incorporating that symmetry
        can significantly improve reconstruction quality and signal recovery.

        The protocol allows the user to define a target symmetry group.
        However, the initial reconstruction process may first proceed in
        asymmetric mode before symmetry is later imposed. This strategy is
        particularly useful because it reduces the risk of enforcing an
        incorrect symmetry too early in the reconstruction process.

        For biological systems with uncertain or ambiguous symmetry, it is
        often safer to begin conservatively and inspect the resulting map
        before applying strong symmetry constraints. Incorrect symmetry
        assignment can artificially distort structural features and may
        hide biologically meaningful asymmetry.

        When symmetry is correctly chosen, however, it can greatly stabilize
        reconstruction and improve map interpretability, especially for
        highly symmetric assemblies such as viral capsids, ring-shaped
        complexes, or oligomeric membrane proteins.

        Classification and Structural Heterogeneity

        The protocol can generate either a single consensus model or
        multiple initial classes. Multiple classes become particularly
        important when the dataset contains conformational variability,
        compositional heterogeneity, or structurally distinct particle
        populations.

        In biological practice, heterogeneous datasets are extremely common.
        Flexible molecular machines, partially assembled complexes, and
        ligand-dependent conformations may all coexist within the same
        experiment. Generating several initial models allows these distinct
        states to begin separating early in the processing workflow.

        However, increasing the number of classes also increases computational
        complexity and may dilute the number of particles contributing to
        each reconstruction. For smaller datasets, excessive classification
        can produce unstable or noisy maps. Users should therefore balance
        biological expectations against the practical limitations of particle
        count and data quality.

        Angular Sampling and Translational Searches

        The protocol allows control over angular and translational search
        parameters during optimization. These settings influence both the
        robustness of reconstruction and the computational cost.

        Fine angular sampling improves orientation precision but requires
        substantially more computation. Coarser sampling is often sufficient
        during early exploratory stages, especially when the objective is
        only to obtain an approximate initial model.

        Translational searches determine how broadly particles may shift
        relative to the evolving reference. Wider ranges are useful when
        particle centering is uncertain, although excessively large search
        spaces may slow convergence and increase the risk of unstable
        solutions.

        In most biological workflows, the default settings provide a good
        balance between robustness and efficiency, particularly for standard
        subtomogram averaging datasets.

        Solvent Flattening and Regularization

        The protocol includes options that stabilize reconstruction by
        controlling noise and solvent behavior. Solvent flattening is
        especially useful for isolated particles embedded within noisy
        tomographic environments because it suppresses irrelevant density
        outside the molecular region.

        Regularization parameters influence how strongly the reconstruction
        is constrained during optimization. Stronger regularization generally
        produces smoother and more stable maps, while weaker regularization
        may preserve additional structural detail at the cost of increased
        sensitivity to noise.

        From a biological perspective, early de novo models should prioritize
        robustness and interpretability rather than maximal detail. Stable
        low-resolution maps provide a much stronger foundation for downstream
        refinement than unstable high-frequency features.

        Computational Considerations

        De novo reconstruction can be computationally demanding because the
        protocol simultaneously estimates particle orientations and the
        evolving reference structure. GPU acceleration substantially improves
        performance and is strongly recommended for practical use with large
        subtomogram datasets.

        The protocol is intended to operate using a single MPI process. This
        limitation reflects the optimization strategy employed during the
        reconstruction procedure and should be considered when planning
        computational resources.

        Particle dimensions are also important. Pseudo-subtomograms should
        have even box sizes to ensure compatibility with the reconstruction
        workflow and avoid downstream processing problems.

        Outputs and Their Interpretation

        After completion, the protocol produces one or more reconstructed
        volumes representing the estimated initial structural states. These
        maps are typically low to intermediate resolution and are intended
        primarily as starting references for subsequent refinement procedures.

        When multiple classes are generated, each output volume represents a
        distinct structural population identified within the dataset. These
        classes may correspond to biologically meaningful conformational
        states, assembly intermediates, or compositional variants.

        The resulting maps should always be inspected visually before further
        processing. Biological plausibility, symmetry consistency, and overall
        structural continuity are important indicators of successful
        reconstruction.

        Practical Recommendations

        In most practical workflows, users should begin with a single class
        and moderate angular sampling to establish whether the dataset can
        produce a stable consensus structure. Once a reliable initial model
        is obtained, additional classes or finer refinement strategies may
        be introduced if structural heterogeneity is suspected.

        Careful particle preparation remains critical. Poorly centered
        particles, insufficient angular coverage, or highly contaminated
        datasets can strongly reduce reconstruction quality regardless of
        optimization settings.

        Symmetry should be applied cautiously, especially during exploratory
        projects involving poorly characterized assemblies. It is often safer
        to verify structural features visually before enforcing strong
        symmetry constraints.

        Final Perspective

        For cryo-electron tomography users, de novo initial model generation
        represents a foundational step in structural interpretation. A robust
        initial reference enables reliable downstream alignment, classification,
        and refinement while minimizing bias introduced by external templates.

        Successful reconstruction depends not only on computational settings
        but also on the biological quality of the dataset itself. Broad angular
        coverage, accurate particle extraction, and realistic expectations of
        achievable resolution are the key factors that determine whether a
        meaningful initial structure can be obtained.
    """

    _label = '3D initial model'
    _possibleOutputs = outputObjects

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    # -------------------------- DEFINE param functions -----------------------

    def _defineParams(self, form):
        self._defineIOParams(form)
        self._defineCTFParams(form)
        self._defineOptimisationParams(form)
        self._defineComputeParams(form)
        self._insertGpuParams(form)
        self._defineAdditionalParams(form)
        form.addParallelSection(threads=0, mpi=1)

    def _defineOptimisationParams(self, form):
        self._insertOptimisationSection(form)
        self._insertVdamMiniBatchesParam(form)
        self._insertRegularisationParam(form)
        self._insertNumOfClassesParam(form)
        self._insertMaskDiameterParam(form)
        self._insertFlattenSolventParam(form)
        helpDeNovo = 'The initial model is always generated in C1 and then aligned to and symmetrized ' \
                     'with the specified point group. If the automatic alignment fails, please manually  ' \
                     'rotate run_itNNN_class001.mrc (NNN is the number of iterations) so that it conforms ' \
                     'the symmetry convention.' + SYMMETRY_HELP_MSG
        self._insertSymmetryParam(form, helpDeNovo)
        self._insertDoInC1AndApplySymLaterParam(form)
        if IS_RELION_50:
            self._insertPriorWidthParam(form)
        self._insertAngularCommonParams(form,
                                        expertLevel=LEVEL_ADVANCED,
                                        angSampling=1,
                                        offsetRange=6,
                                        offsetStep=2)

    # -------------------------- INSERT steps functions -----------------------
    def _insertAllSteps(self):
        self._insertFunctionStep(self.convertInputStep, needsGPU=False)
        self._insertFunctionStep(self.generateDeNovo3DModel, needsGPU=True)
        self._insertFunctionStep(self.alignSymmetry, needsGPU=False)
        self._insertFunctionStep(self.createOutputStep, needsGPU=False)

    # -------------------------- STEPS functions ------------------------------
    def generateDeNovo3DModel(self):
        nMpi = self.numberOfMpi.get()
        Plugin.runRelionTomo(self,
                             getProgram('relion_refine', nMpi),
                             self._genInitModelCommand(),
                             numberOfMpi=nMpi)

    def alignSymmetry(self):
        if self.numberOfClasses.get() == 1:
            Plugin.runRelionTomo(self,
                                 'relion_align_symmetry',
                                 self._genApplySymCmd())
        else:
            for i in range(self.numberOfClasses.get()):
                Plugin.runRelionTomo(self,
                                     'relion_align_symmetry',
                                     self._genApplySymCmd(classIndex=i + 1))

    def createOutputStep(self):
        inRelionParticlesPointer = self.getInputParticles(returnPointer=True)
        if self.numberOfClasses.get() == 1:
            vol = self._createOutputModel(inRelionParticlesPointer)
            self._defineOutputs(**{outputObjects.average.name: vol})
            self._defineSourceRelation(inRelionParticlesPointer, vol)
        else:
            avgSet = SetOfVolumes.create(self._getPath(), template='averages%s.sqlite')
            avgSet.copyInfo(self.getInputParticles())
            for i in range(self.numberOfClasses.get()):
                vol = self._createOutputModel(inRelionParticlesPointer, classIndex=i + 1)
                avgSet.append(vol)

            self._defineOutputs(**{outputObjects.averages.name: avgSet})
            self._defineSourceRelation(inRelionParticlesPointer, avgSet)

    # -------------------------- INFO functions -------------------------------
    def _validate(self):
        errorMsg = []
        nMpi = self.numberOfMpi.get()
        if nMpi > 1:
            errorMsg.append('The initial volume can only run using 1 MPI.')

        if self.inReParticles.get().getBoxSize() % 2 != 0:
            errorMsg.append('The dimensions of the extracted pseudo-subtomograms '
                            'must be even!')
        return errorMsg

    # --------------------------- UTILS functions -----------------------------
    def _genInitModelCommand(self) -> str:
        # Common parameters from base protocol
        cmd = self._genBaseCommand(useOptimizationSet=False)

        # Initial model specific commands
        cmd += ' --denovo_3dref --grad --zero_mask --auto_sampling --pad 1 '
        #   Optimisation args
        cmd += '--iter %i ' % self.nVdamMiniBatches.get()
        cmd += '--tau2_fudge %d ' % self.regularisation.get()
        cmd += '--K %i ' % self.numberOfClasses.get()
        if self.flattenSolvent.get():
            cmd += '--flatten_solvent '
        if self.doInC1AndApplySymLater.get():
            cmd += '--sym C1 '
        else:
            cmd += '--sym %s ' % self.symmetry.get()
        cmd += '--healpix_order %i ' % self.angularSamplingDeg.get()
        cmd += '--offset_step %i ' % self.offsetSearchStepPix.get()
        cmd += '--offset_range %d ' % self.offsetSearchRangePix.get()
        if IS_RELION_50:
            cmd += '--sigma_tilt %i ' % self.priorWidthTiltAngle.get()
        return cmd

    def _genApplySymCmd(self, classIndex: Union[int, None] = None) -> str:
        cmd = '--o %s ' % self._getInitialModelOutFn(classIndex=classIndex)
        classIndex = 1 if classIndex is None else classIndex
        cmd += '--i %s ' % self._getExtraPath(self._getModelName(classIndex))
        if self.doInC1AndApplySymLater.get() and 'c1' not in self.symmetry.get().lower():
            cmd += '--sym %s ' % self.symmetry.get()
        else:
            cmd += '--sym C1 '
        cmd += '--apply_sym --select_largest_class '
        return cmd

    def _getModelName(self, classIndex: Union[int, None]) -> str:
        """generate the name of the volume following this pattern _it002_model.star"""
        classIndexStr = f'{classIndex:03d}' if classIndex is not None else ''
        return f'_it{self.nVdamMiniBatches.get():03d}_class{classIndexStr}.mrc'

    def _getInitialModelOutFn(self, classIndex: Union[int, None] = None) -> str:
        classIndexStr = f'_{classIndex:03d}' if classIndex is not None else ''
        return self._getExtraPath(f'initial_model{classIndexStr}.mrc')

    def _createOutputModel(self,
                           inRelionParticlesPointer: Pointer,
                           classIndex: Union[int, None] = None) -> Volume:
        vol = Volume()
        iniModelFile = self._getInitialModelOutFn(classIndex=classIndex)
        fixVolume(iniModelFile)  # Fix header to make it interpreted as volume instead of a stack by xmipp
        vol.setFileName(iniModelFile)
        vol.setSamplingRate(inRelionParticlesPointer.get().getSamplingRate())
        return vol
