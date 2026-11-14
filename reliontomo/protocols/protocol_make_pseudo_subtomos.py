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
from pyworkflow.protocol import FloatParam, BooleanParam, GE
from pyworkflow.protocol.constants import LEVEL_ADVANCED
from reliontomo import Plugin
from reliontomo.objects import RelionSetOfPseudoSubtomograms
from reliontomo.protocols.protocol_base_make_pseusosubtomos_and_rec_particle import \
    ProtRelionMakePseudoSubtomoAndRecParticleBase
from reliontomo.utils import getProgram
from tomo.protocols import ProtTomoBase


class outputObjects(Enum):
    relionParticles = RelionSetOfPseudoSubtomograms


class ProtRelionMakePseudoSubtomograms(ProtRelionMakePseudoSubtomoAndRecParticleBase, ProtTomoBase):
    """
    Make pseudo-subtomograms

    Pseudo-subtomograms do not aim to accurately represent the scattering potential of the underlying particles.
    Instead, they serve as a practical means to implement an approximation to the 2D approach within
    the existing RELION framework. In the original RELION4 article an accurate definition is given.
    More info: https://doi.org/10.7554/eLife.83724

    AI Generated:

    Make Pseudo-Subtomograms (ProtRelionMakePseudoSubtomograms) - User Manual
        Overview

        The Make Pseudo-Subtomograms protocol generates pseudo-subtomograms
        from tomographic particle data for use within RELION subtomogram
        averaging workflows. Rather than reconstructing fully physical 3D
        particle densities, pseudo-subtomograms provide an efficient
        approximation that allows tomography data to be processed using
        RELION refinement and classification strategies originally designed
        for single-particle analysis.

        In biological cryo-electron tomography projects, this protocol acts
        as a bridge between tilt-series information and downstream alignment,
        refinement, and classification procedures. It converts particle
        observations from multiple tilt images into a form that can be
        interpreted consistently within RELION iterative optimization
        workflows while preserving the angular and contrast transfer
        information required for high-resolution subtomogram averaging.

        Biological Purpose and Context

        Pseudo-subtomograms are commonly used when studying macromolecular
        complexes directly inside cells, organelles, membranes, or native
        biological environments. In these workflows, particles are observed
        across multiple tilted projections rather than through isolated
        single-particle micrographs. This protocol organizes those
        observations into particle-centered volumetric representations that
        support alignment and averaging while remaining computationally
        practical for large tomography datasets.

        From a biological perspective, the protocol is especially useful for
        studying heterogeneous molecular assemblies, membrane-associated
        proteins, ribosomes in situ, viral particles, or cellular complexes
        embedded within crowded environments. By preparing particles for
        RELION tomography refinement, the protocol contributes directly to
        improving structural interpretability and averaging consistency.

        Inputs and Reconstruction Workflow

        The protocol requires reconstructed particle information and the
        associated tomographic metadata describing how each particle was
        observed throughout the tilt series. These inputs are used to build
        pseudo-subtomograms that integrate contrast transfer information and
        sampling statistics from the original tilt images.

        The resulting pseudo-subtomograms are not intended to represent
        final biological densities for visualization or interpretation.
        Instead, they are optimized for iterative refinement and statistical
        comparison inside RELION. Their primary role is to support accurate
        alignment and classification during downstream processing.

        Cone Weighting and Membrane Suppression

        An important optional feature of this protocol is cone weighting,
        which suppresses Fourier-space information along a selected direction.
        This functionality is particularly relevant for membrane-associated
        particles. In many tomographic datasets, membrane signal can dominate
        the alignment process because planar membranes generate strong and
        highly directional Fourier contributions.

        By reducing the influence of this directional membrane information,
        the protocol allows alignment to focus more strongly on the embedded
        macromolecular complex itself rather than on the surrounding lipid
        bilayer. This can significantly improve refinement stability for
        membrane proteins, membrane-bound assemblies, and complexes attached
        to curved cellular surfaces.

        The cone angle determines how broadly the suppression region extends.
        Biologically, wider cone angles may be appropriate when membrane
        orientations vary substantially or when particle orientations remain
        uncertain. Narrower angles preserve more information but may provide
        less effective suppression of dominant membrane features.

        Precision and Storage Considerations

        The protocol optionally allows reduced-precision storage formats in
        order to decrease disk usage. This can be particularly valuable in
        large cryo-electron tomography projects where thousands or millions
        of pseudo-subtomograms may be generated. Although lower-precision
        storage is often sufficient for iterative refinement workflows,
        users should verify compatibility with downstream visualization and
        processing software before adopting compressed representations in
        production pipelines.

        Outputs and Their Interpretation

        The protocol produces a set of pseudo-subtomograms prepared for
        RELION tomography refinement and classification workflows. These
        outputs retain the particle relationships and metadata needed for
        downstream subtomogram averaging while providing an efficient
        intermediate representation for computational processing.

        Biologically, the generated pseudo-subtomograms should be viewed as
        analytical constructs designed for refinement rather than direct
        physical reconstructions. Their value lies in enabling accurate
        particle alignment, classification, and averaging across complex
        tomographic datasets.

        Practical Recommendations

        In routine cryo-electron tomography workflows, users commonly apply
        this protocol after particle picking and initial coordinate
        refinement. Careful attention should be given to the quality of the
        input particle orientations and tomographic metadata, since errors at
        this stage can propagate into downstream averaging results.

        For membrane-associated systems, cone weighting often provides a
        major improvement in alignment robustness by reducing bias from the
        surrounding bilayer. However, excessive suppression may also remove
        biologically meaningful signal, so parameter selection should be
        guided by both structural expectations and visual inspection of
        refinement outcomes.

        Final Perspective

        Pseudo-subtomogram generation is a central preparatory step in
        RELION tomography processing because it transforms raw tomographic
        particle observations into a form suitable for iterative refinement
        and classification. By balancing computational efficiency with the
        preservation of biologically relevant signal, this protocol enables
        high-resolution subtomogram averaging workflows for complex cellular
        and in situ structural biology applications.
    """

    _label = 'Make pseudo-subtomograms'
    _possibleOutputs = outputObjects

    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        super()._defineParams(form)
        form.addSection(label='Reconstruct pseudo-Subtomograms')
        self._defineCommonRecParams(form)
        form.addParam('applyConeWeight', BooleanParam,
                      label='Apply cone weight?',
                      default=False,
                      help='Down weight a cone in Fourier space along the Z axis (as defined by the coordinate system '
                           'of the particle). This is useful for particles embedded in a membrane, as it can prevent '
                           'the alignment from being driven by the membrane signal (the signal of a planar membrane is '
                           'localised within one line in 3D Fourier space). Note that the coordinate system of a '
                           'particle is given by both the subtomogram orientation (if defined) and the particle '
                           'orientation. This allows the user to first obtain a membrane-driven alignment, and to then '
                           'specifically suppress the signal in that direction.')
        form.addParam('coneAngle', FloatParam,
                      label='Cone angle (deg.)',
                      condition='applyConeWeight',
                      validator=[GE(0)],
                      default=10,
                      help='It is the (full) opening angle of the cone to be suppressed, given in degrees. This angle '
                           'should  include both the uncertainty about the membrane orientation and its variation '
                           'across the region represented in the subtomogram.')
        form.addParam('outputInFloat16', BooleanParam,
                      label='Write output in float16?',
                      expertLevel=LEVEL_ADVANCED,
                      default=False,
                      help='If set to Yes, this program will write output images in float16 MRC format. This will '
                           'save a factor of two in disk space compared to the default of writing in float32. Note '
                           'that RELION and CCPEM will read float16 images, but other programs may not (yet) do so.')
        self._defineExtraParams(form)

    # -------------------------- INSERT steps functions -----------------------
    def _insertAllSteps(self):
        self._insertFunctionStep(self.convertInputStep, needsGPU=False)
        self._insertFunctionStep(self.relionMakePseudoSubtomos, needsGPU=False)
        self._insertFunctionStep(self.createOutputStep, needsGPU=False)

    # -------------------------- STEPS functions ------------------------------
    def relionMakePseudoSubtomos(self):
        Plugin.runRelionTomo(self, getProgram('relion_tomo_subtomo', self.numberOfMpi.get()),
                             self._genMakePseudoSubtomoCmd(), numberOfMpi=self.numberOfMpi.get())

    def createOutputStep(self):
        psubtomoSet = super().createOutputStep()
        self._defineOutputs(**{outputObjects.relionParticles.name: psubtomoSet})
        self._defineSourceRelation(self.inReParticles.get(), psubtomoSet)

    # -------------------------- INFO functions -------------------------------

    # # --------------------------- UTILS functions -----------------------------
    def _genMakePseudoSubtomoCmd(self):
        cmd = self._genCommonCmd()
        cmd += '--o %s ' % self._getExtraPath()
        if self.applyConeWeight.get():
            cmd += '--cone_weight '
            cmd += '--cone_angle %.2f ' % self.coneAngle.get()
        if self.outputInFloat16.get():
            cmd += '--float16 '
        return cmd
