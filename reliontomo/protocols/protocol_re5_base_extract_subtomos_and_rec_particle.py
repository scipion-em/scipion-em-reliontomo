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
from pyworkflow.utils import createLink
from reliontomo.constants import IN_PARTICLES_STAR, IN_TOMOS_STAR
from reliontomo.protocols.protocol_base_relion import ProtRelionTomoBase, IS_RELION_50


class ProtRelion5ExtractSubtomoAndRecParticleBase(ProtRelionTomoBase):
    """
    Provides the common foundation for reconstructing particle-centered
    tomographic regions and generating pseudo-subtomograms from cryo-electron
    tomography data. The protocol is intended for workflows in which particles
    identified within tomograms must be re-extracted into localized 3D volumes
    suitable for subtomogram refinement, classification, structural inspection,
    or downstream Relion-based analysis.

    AI Generated:

    Reconstruct Particle and Generate Pseudo-Subtomograms Base Protocol
    (ProtRelion5ExtractSubtomoAndRecParticleBase) - User Manual

        Overview

        This protocol serves as the shared reconstruction framework for
        particle-centered tomographic extraction workflows in Relion 5. Its
        main objective is to reconstruct localized three-dimensional regions
        around particles identified within tomograms and to prepare these
        regions as pseudo-subtomograms that can be used in subtomogram
        averaging and refinement pipelines.

        In cryo-electron tomography, particles are often first identified as
        coordinates or as refined particle trajectories derived from previous
        alignment procedures. Once these particle positions are known, users
        typically need to reconstruct consistent local volumes around each
        particle in order to improve alignment precision, perform classification,
        or inspect structural heterogeneity. This protocol provides the common
        reconstruction logic required for those tasks.

        Biological Context and Typical Applications

        For biological users, localized reconstruction around particles is an
        essential step in high-resolution subtomogram analysis. Rather than
        analyzing the entire tomogram, which may contain large amounts of noise,
        contamination, or unrelated structures, the workflow isolates the
        molecular region surrounding each particle and reconstructs it into a
        compact 3D representation.

        This approach is widely used for membrane proteins embedded in native
        cellular environments, ribosomes inside cells, viral assemblies,
        cytoskeletal complexes, and many other macromolecular systems studied
        directly in situ. By reconstructing only the relevant particle-centered
        region, downstream alignment and averaging become more computationally
        efficient and biologically interpretable.

        Inputs and Reconstruction Workflow

        The protocol operates on particle information associated with tomograms.
        These particles may originate from coordinate picking, previous
        subtomogram refinement, or trajectory-based motion analysis. The
        workflow assumes that the tomographic metadata and particle metadata
        remain properly associated so that reconstruction occurs in the correct
        spatial context.

        During execution, the protocol prepares the information required for
        Relion-based localized reconstruction. This includes particle metadata,
        tomogram references, and optional trajectory information describing
        particle motion or refinement history. The resulting reconstruction is
        then generated around each particle position according to the selected
        extraction parameters.

        Binning and Downsampling

        One of the most important reconstruction parameters is the binning
        factor. Binning reduces the effective sampling of the tilt-series images
        before reconstruction, which decreases computational cost and storage
        requirements. Larger binning factors accelerate processing and are often
        useful during exploratory analyses or initial refinement stages.

        However, biological users should recognize that excessive downsampling
        may remove high-resolution structural information. Small protein
        complexes, fine membrane features, or flexible domains may become poorly
        resolved if the sampling becomes too coarse. For this reason, lower
        binning values are generally preferred for final refinement and
        high-resolution interpretation.

        Reconstruction Box Size

        The reconstruction box size determines the physical region surrounding
        each particle that will be reconstructed into a three-dimensional map.
        Choosing an appropriate box size is biologically important because it
        controls how much contextual information is preserved around the target
        structure.

        Small box sizes improve computational efficiency and are suitable for
        compact particles with limited surrounding density. Larger box sizes are
        advantageous when the particle interacts with membranes, neighboring
        complexes, or flexible assemblies that may contribute biologically
        meaningful signal.

        Users should avoid selecting boxes that are too small, since important
        structural information or delocalized high-frequency signal may be
        excluded. Conversely, excessively large boxes increase storage and
        processing demands without necessarily improving interpretability.

        Cropped Box Sizes and Refinement Efficiency

        The protocol also supports generating cropped pseudo-subtomograms after
        reconstruction. This allows the reconstruction step to preserve a large
        contextual region while producing smaller output volumes optimized for
        refinement.

        In practical cryo-ET workflows, this strategy is highly valuable because
        refinement and classification of large subtomograms can become extremely
        computationally expensive. Cropping reduces memory usage and accelerates
        downstream alignment procedures while still preserving the most relevant
        structural region.

        Particle Trajectories and Motion Information

        In advanced workflows, particles may contain trajectory information
        derived from motion refinement or particle tracking procedures. These
        trajectories help maintain accurate particle positioning throughout the
        reconstruction process and can improve the biological consistency of the
        reconstructed subtomograms.

        This becomes particularly important for in situ datasets where beam
        induced motion, local deformation, or sample instability may affect
        particle positioning. Incorporating trajectory information can therefore
        improve alignment quality and contribute to higher-resolution averages.

        Outputs and Downstream Applications

        The resulting pseudo-subtomograms are intended for use in Relion 5
        subtomogram refinement and classification workflows. These outputs
        provide localized 3D particle representations that preserve the spatial
        orientation and reconstruction context required for accurate alignment.

        Biologically, the generated subtomograms may be used to identify
        conformational states, improve structural resolution, analyze molecular
        variability, or study macromolecular organization directly inside cells
        and tissues.

        Practical Recommendations

        For exploratory processing, users often begin with moderate binning and
        relatively small cropped box sizes in order to accelerate testing and
        parameter optimization. Once particle quality and alignment stability
        have been confirmed, reconstruction can be repeated using lower binning
        and larger contextual boxes for higher-quality structural analysis.

        When studying membrane-associated complexes or assemblies embedded in
        crowded cellular environments, it is generally beneficial to preserve
        enough surrounding density to maintain meaningful biological context.
        However, refinement-focused workflows usually benefit from tighter
        cropped boxes centered on the particle of interest.

        Final Perspective

        Localized particle reconstruction is one of the central operations in
        subtomogram averaging workflows because it bridges tomographic particle
        detection with high-resolution structural refinement. Careful selection
        of binning, reconstruction size, and contextual coverage strongly
        influences both computational efficiency and biological interpretability.

        By providing a unified reconstruction foundation for pseudo-subtomogram
        generation, this protocol supports scalable and biologically meaningful
        cryo-electron tomography analysis within Relion 5 environments.
    """

    _label = None

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.isInSetOf3dCoords = None

    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        self._defineCommonInputParams(form)
        self._insertBinThreadsParam(form)
        form.addParallelSection(threads=0, mpi=3)

    @staticmethod
    def _defineCommonRecParams(form):
        form.addParam('binningFactor', FloatParam,
                      label='Binning factor (downsampling)',
                      default=1,
                      validator=GE(0),
                      allowsNull=False,
                      important=True,
                      help='The tilt series images will be binned by this (real-valued) factor and then '
                           'reconstructed in the specified box size above. Note that thereby the '
                           'reconstructed region becomes larger when specifying binning factors larger than one. '
                           'This does not alter the box size.')
        form.addParam('boxSize', IntParam,
                      label='Box size (px)',
                      validator=GE(0),
                      allowsNull=False,
                      important=True,
                      help='Box size, in pixels, of the reconstruction. Note that this is independent of the '
                           'box size used to refine the particle. This allows the user to construct a 3D map of '
                           'arbitrary size to gain an overview of the structure surrounding the particle. A '
                           'sufficiently large box size also allows more of the high-frequency signal to be '
                           'captured that has been delocalized by the CTF.')
        form.addParam('croppedBoxSize', IntParam,
                      label="Cropped box size (px)",
                      allowsNull=True,
                      important=True,
                      help='The resulting pseudo subtomograms are cropped to this size. A smaller box size '
                           'allows the (generally expensive) refinement using relion_refine to proceed more rapidly.')

    # -------------------------- STEPS functions ------------------------------
    def convertInputStep(self):
        inParticles = self.getInputParticles()
        # Generate the file particles.star
        self.genInStarFile(are2dParticles=inParticles.are2dStacks())
        # Link the file tomograms.star
        # The tomograms file will exist and be stored as an attribute of the set, having been updated if a new one is
        # generated, like in the protocol bayesian polishing
        createLink(inParticles.getTomogramsStar(), self._getExtraPath(IN_TOMOS_STAR))
        # Tilt-series star files:
        # The tilt-series star files will exist and their corresponding path will be provided by the file tomograms.star

    # --------------------------- UTILS functions -----------------------------
    @classmethod
    def isDisabled(cls):
        """ Return True if this Protocol is disabled.
        Disabled protocols will not be offered in the available protocols."""
        return False if IS_RELION_50 else True

    def _genCommonExtractAndRecCmd(self):
        cmd = [f'--p {self._getExtraPath(IN_PARTICLES_STAR)}',
               f'--t {self._getExtraPath(IN_TOMOS_STAR)}',
               f'--o {self._getExtraPath()}',
               f"--b {self.boxSize.get()}",
               f"--crop {self.croppedBoxSize.get()}",
               f"--bin {self.binningFactor.get():.1f}",
               f"--j {self.binThreads.get()}",
               self._genExtraParamsCmd()]
        inParticles = self.getInputParticles()
        if not self.isInputSetOf3dCoords():
            trajectoriesFile = inParticles.getTrajectoriesStar()
            if trajectoriesFile:
                cmd.append(f'--mot {trajectoriesFile}')
        return ' '.join(cmd)
