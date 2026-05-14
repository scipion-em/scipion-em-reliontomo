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
from reliontomo.protocols.protocol_base_relion import ProtRelionTomoBase, IS_RELION_50


class ProtRelionMakePseudoSubtomoAndRecParticleBase(ProtRelionTomoBase):
    """
    Reconstruct particle and make pseudo-subtomograms base class

    AI Generated:

    Make Pseudo-Subtomograms and Reconstruct Particle Base
    (ProtRelionMakePseudoSubtomoAndRecParticleBase) - User Manual

        Overview

        This protocol provides the foundational framework for generating
        pseudo-subtomograms and reconstructed particle volumes within
        cryo-electron tomography workflows. Its primary objective is to
        transform aligned tomographic particle information into
        reconstructed volumetric representations that can later be used
        for subtomogram averaging, refinement, classification, and
        high-resolution structural interpretation.

        In practical biological workflows, pseudo-subtomograms serve as
        intermediate representations that allow tomography data to be
        processed using refinement strategies originally developed for
        single-particle cryo-electron microscopy. This approach enables
        the integration of tomographic information into highly optimized
        refinement pipelines while preserving the geometrical and
        experimental relationships associated with tilt-series
        acquisition.

        Biological Purpose

        Cryo-electron tomography datasets are intrinsically complex
        because each particle is reconstructed from multiple tilted
        projection images acquired under varying imaging conditions.
        This protocol establishes the infrastructure required to convert
        those experimental observations into localized 3D particle
        volumes suitable for downstream refinement and analysis.

        From a biological perspective, this reconstruction stage is
        critical because it defines the quality and interpretability of
        all subsequent analyses. Accurate pseudo-subtomograms allow
        users to identify structural details, classify conformational
        states, and refine macromolecular assemblies embedded within
        native cellular environments.

        Inputs and Reconstruction Workflow

        The protocol requires pseudo-subtomogram metadata together with
        the associated tomographic information describing the original
        tilt-series geometry and particle relationships. These inputs
        are assumed to originate from earlier tomographic processing
        stages where particles have already been detected, extracted,
        and aligned.

        During reconstruction, the protocol generates volumetric
        particle representations using configurable spatial sampling and
        reconstruction dimensions. The resulting pseudo-subtomograms are
        designed to preserve biologically meaningful signal while
        remaining computationally manageable for downstream refinement.

        The workflow is intended to integrate naturally into iterative
        refinement pipelines, where reconstructed particles are
        repeatedly improved through alignment, classification, and local
        correction procedures.

        Reconstruction Box Size

        One of the most biologically important parameters is the
        reconstruction box size. This parameter determines the spatial
        region surrounding the particle that will be reconstructed into
        the volumetric representation.

        Larger box sizes allow the inclusion of surrounding structural
        context, which can be valuable for membrane proteins, large
        macromolecular assemblies, or complexes interacting with nearby
        cellular components. A sufficiently large reconstruction volume
        also helps preserve high-frequency information that may become
        spatially delocalized by the contrast transfer function.

        However, excessively large reconstruction boxes substantially
        increase memory usage, computational demand, and storage
        requirements. In routine workflows, the selected box size should
        balance structural completeness with computational efficiency.

        Cropped Reconstruction Volumes

        The protocol also supports generation of cropped particle
        volumes for downstream refinement. Cropping allows the user to
        reconstruct a larger contextual region while later retaining a
        smaller computationally efficient region for intensive
        refinement procedures.

        This strategy is particularly useful in high-resolution
        tomography workflows, where large contextual reconstructions may
        improve signal recovery but smaller particle boxes are preferred
        during iterative refinement. Care must be taken to ensure that
        the cropped region still fully contains the biologically
        relevant structure.

        Binning and Sampling Considerations

        The protocol allows reconstruction at different binning levels,
        which directly affects voxel sampling and computational cost.
        Binning reduces the effective image resolution while improving
        processing speed and reducing memory consumption.

        In exploratory analyses or early refinement stages, moderate
        binning is often advantageous because it accelerates processing
        and stabilizes alignment. For final high-resolution analyses,
        however, lower binning or full-resolution reconstructions are
        generally preferred.

        Importantly, increasing the binning factor changes the spatial
        sampling of the reconstructed region without necessarily
        reducing the overall reconstruction dimensions. This allows
        broader structural context to remain visible while maintaining
        manageable computational requirements.

        Outputs and Their Interpretation

        The protocol produces reconstructed pseudo-subtomograms together
        with updated metadata suitable for subsequent Relion tomography
        workflows. These outputs can be directly used for refinement,
        classification, motion correction, or local optimization
        procedures.

        Biologically, the resulting reconstructed particles provide the
        structural foundation for all downstream analyses. The quality
        of these reconstructions strongly influences the achievable
        resolution, interpretability of flexible regions, and reliability
        of conformational classification.

        Practical Recommendations

        In most workflows, users should initially select reconstruction
        parameters that prioritize stability and computational
        efficiency. Moderate binning and carefully chosen cropped box
        sizes often provide the best compromise during early iterative
        refinement.

        Larger reconstruction boxes are beneficial when surrounding
        structural context is biologically relevant or when strong CTF
        delocalization effects are expected. However, cropped volumes
        should remain sufficiently large to fully contain the particle
        throughout all refinement stages.

        Users should avoid selecting identical reconstruction and
        cropped box sizes when later refinement steps expect reduced
        particle dimensions. Consistent parameter selection throughout
        the tomography workflow helps prevent downstream incompatibility
        and improves refinement robustness.

        Final Perspective

        Reconstruction of pseudo-subtomograms is one of the central
        stages in cryo-electron tomography processing because it
        converts experimental tilt-series information into particle-
        centered volumetric representations suitable for structural
        analysis. By balancing reconstruction size, cropping strategy,
        and sampling resolution, this protocol provides the essential
        foundation for accurate and biologically meaningful tomographic
        refinement workflows.
    """

    _label = None

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        self._defineCommonInputParams(form)
        self._insertBinThreadsParam(form)
        form.addParallelSection(threads=0, mpi=3)

    @staticmethod
    def _defineCommonRecParams(form):
        form.addParam('boxSize', IntParam,
                      label='Box size (px)',
                      validator=GE(0),
                      important=True,
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
                           ' allows the (generally expensive) refinement using relion_refine to proceed more rapidly.')
        form.addParam('binningFactor', FloatParam,
                      label='Binning factor (downsampling)',
                      default=1,
                      validator=GE(0),
                      allowsNull=False,
                      help='The tilt series images will be binned by this (real-valued) factor and then '
                           ' reconstructed in the specified box size above. Note that thereby the '
                           ' reconstructed region becomes larger when specifying binning factors larger than one.'
                           ' This does not alter the box size.')

    # -------------------------- INSERT steps functions -----------------------
    def _insertAllSteps(self):
        pass

    def convertInputStep(self):
        self.genInStarFile()

    def createOutputStep(self):
        return self.genRelionParticles(binningFactor=self.binningFactor.get(), boxSize=self.croppedBoxSize.get())

    # -------------------------- INFO functions -------------------------------
    def _warnings(self):
        warnMsg = []
        boxSize = self.boxSize.get()
        croppedBoxSize = self.croppedBoxSize.get()
        if boxSize == croppedBoxSize:
            warnMsg.append('Setting the same value to the Box size and the Cropped box size may cause errors in '
                           'later steps of the refinement cycle.')
        elif boxSize < croppedBoxSize:
            warnMsg.append(f"The Box size [{boxSize} px] should be lower than the Cropped box size [{croppedBoxSize} "
                           f"px]. Please check these parameters' help to get more detailed information.")
        return warnMsg

    # --------------------------- UTILS functions -----------------------------
    @classmethod
    def isDisabled(cls):
        """ Return True if this Protocol is disabled.
        Disabled protocols will not be offered in the available protocols."""
        return True if IS_RELION_50 else False

    def _genCommonCmd(self):
        inRelionParticles = self.getInputParticles()
        cmd = ''

        # Cancel this for now
        self.info("Using optimization_set: %s" % inRelionParticles.filesMaster)
        cmd += '--i %s ' % inRelionParticles.filesMaster

        # This would be either the particles start file in the set or a new generated one from the subtomo set.
        cmd += '--p %s ' % self.getOutStarFileName()

        if inRelionParticles.getTrajectoriesStar():
            cmd += '--mot %s ' % inRelionParticles.getTrajectoriesStar()

        cmd += '--b %i ' % self.boxSize.get()
        cmd += '--crop %i ' % self.croppedBoxSize.get()
        cmd += '--bin %.1f ' % self.binningFactor.get()
        cmd += '--j %i ' % self.binThreads.get()
        cmd += self._genExtraParamsCmd()
        return cmd
