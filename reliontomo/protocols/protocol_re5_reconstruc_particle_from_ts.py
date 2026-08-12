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
from enum import Enum
from os.path import exists

from pwem.convert.headers import fixVolume
from pyworkflow.protocol import StringParam, FloatParam
from reliontomo import Plugin
from reliontomo.constants import SYMMETRY_HELP_MSG
from reliontomo.protocols.protocol_re5_base_extract_subtomos_and_rec_particle import \
    ProtRelion5ExtractSubtomoAndRecParticleBase
from tomo.objects import AverageSubTomogram


class outputObjects(Enum):
    average = AverageSubTomogram


class ProtRelion5ReconstructParticle(ProtRelion5ExtractSubtomoAndRecParticleBase):
    """
    Reconstructs and averages subtomographic particle information extracted
    from tilt series data in RELION 5. The protocol generates a 3D map
    representing the consensus structure of the analyzed particles and is
    designed for cryo-electron tomography workflows where aligned particle
    projections are combined into a refined volumetric reconstruction.

    AI Generated:

    Reconstruct Particle Relion 5 (ProtRelion5ReconstructParticle) - User Manual
        Overview

        The Reconstruct Particle Relion 5 protocol performs subtomogram
        reconstruction and averaging using RELION 5 tomography tools. Its
        main purpose is to combine particle information extracted from tilt
        series experiments into a single three-dimensional density map that
        represents the average structure of the analyzed specimen.

        In cryo-electron tomography workflows, this protocol is commonly used
        after particle extraction and alignment steps have already determined
        the orientations and positions of individual particles. The resulting
        reconstruction provides an interpretable 3D map that can be used for
        structural analysis, visualization, classification, or subsequent
        refinement procedures.

        Biological Context

        Subtomogram averaging is particularly useful for studying molecular
        complexes in their native cellular or in situ environments. Unlike
        single-particle analysis, cryo-electron tomography preserves spatial
        context and allows biological assemblies to be analyzed directly
        inside cells, organelles, membranes, or pleomorphic structures.

        By averaging many aligned particles together, the protocol improves
        the signal-to-noise ratio and enhances structural features that may
        not be visible in individual tomograms. This is especially important
        in tomography datasets, where acquisition limitations and electron
        dose constraints typically produce noisier data than conventional
        single-particle experiments.

        Inputs and Reconstruction Workflow

        The protocol expects particles that already contain orientation and
        alignment information derived from previous tomography processing
        steps. These particles are combined into a consensus reconstruction
        that represents the average structural state of the selected dataset.

        The reconstruction process operates directly on tomography-derived
        particle information and supports parallel execution for efficient
        processing of large datasets. This makes the protocol suitable for
        both exploratory analyses and high-throughput cryo-ET studies.

        Symmetry Considerations

        Users may specify a symmetry group for the reconstruction. Symmetry
        can significantly improve the quality of the final map when the
        biological assembly possesses known rotational or point-group
        symmetry. Applying symmetry increases the effective amount of signal
        contributing to the reconstruction and often improves resolution and
        map interpretability.

        However, symmetry should only be imposed when biologically justified.
        Incorrect symmetry assignment may distort structural features or mask
        meaningful asymmetry. For complexes with unknown or flexible
        organization, using C1 symmetry is generally the safest starting
        point.

        Wiener Filtering and Signal Preservation

        The protocol optionally supports Wiener filtering through the use of
        a signal-to-noise ratio estimate. This filtering can help stabilize
        reconstructions and reduce amplification of noise in poorly sampled
        frequency regions.

        From a biological perspective, conservative filtering strategies are
        usually preferable. Excessively aggressive filtering may suppress
        high-resolution information and obscure subtle structural features
        that are important for interpretation or downstream refinement.
        Conversely, moderate filtering can improve the visual clarity and
        robustness of the resulting density map.

        Parallelization and Computational Considerations

        Tomographic reconstruction is computationally demanding due to the
        large size of tomograms and the complexity of subtomogram averaging.
        The protocol is designed to operate in distributed computing
        environments using MPI-based execution, allowing efficient handling
        of large particle datasets.

        Biological users working with extensive in situ datasets should
        ensure that sufficient computational resources are available,
        particularly memory resources. Large reconstructions may require
        substantial memory allocation during averaging and Fourier-space
        processing.

        Outputs and Interpretation

        The primary output is an averaged subtomogram volume representing the
        consensus structure of the reconstructed particles. When available,
        half maps are also preserved, enabling downstream validation and
        resolution estimation procedures.

        The reconstructed map can be used for visualization, segmentation,
        docking of atomic models, structural interpretation, or additional
        refinement workflows. Because the reconstruction represents an
        average over many particles, the resulting density emphasizes common
        structural features while reducing random noise contributions.

        Biological interpretation should nevertheless consider the potential
        presence of structural heterogeneity. Flexible regions, compositional
        variability, or mixed conformational states may appear blurred or
        weakened in the averaged reconstruction.

        Practical Recommendations

        In most biological workflows, reconstruction quality strongly depends
        on the accuracy of the upstream alignment and particle extraction
        procedures. Before interpreting structural details, users should
        visually inspect the resulting map and verify that major features are
        biologically plausible and consistent with prior knowledge.

        Applying symmetry should be approached cautiously, especially for
        complexes embedded in membranes or cellular environments where
        asymmetry may be functionally important. Similarly, Wiener filtering
        parameters should be selected conservatively to avoid suppressing
        relevant structural information.

        For large cryo-ET datasets, distributed execution with multiple MPI
        processes is generally recommended to ensure efficient reconstruction
        performance and manageable processing times.

        Final Perspective

        The Reconstruct Particle Relion 5 protocol provides a practical and
        biologically meaningful framework for generating consensus 3D maps
        from subtomogram datasets. By combining aligned particles into a
        unified reconstruction, it enables structural interpretation of
        macromolecular assemblies directly within their native tomographic
        environments.

        Careful consideration of symmetry, filtering, and data quality is
        essential for obtaining biologically reliable reconstructions that
        accurately reflect the structural organization of the studied system.
    """

    _label = 'Reconstruct particle Relion 5'
    _possibleOutputs = outputObjects

    def __init__(self, **args):
        super().__init__(**args)

    # -------------------------- DEFINE param functions -----------------------

    def _defineParams(self, form):
        super()._defineParams(form)
        form.addSection(label='Average')
        self._defineCommonRecParams(form)
        form.addParam('symmetry', StringParam,
                      label='Symmetry group',
                      default='C1',
                      help=SYMMETRY_HELP_MSG)
        form.addParam('snrWiener', FloatParam,
                      label='Apply a Wiener filter with this SNR',
                      default=0,
                      help='If set to a positive value, apply a Wiener filter with this signal-to-noise ratio. If '
                           'omitted, the reconstruction will use a heuristic to prevent divisions by excessively '
                           'small numbers. Please note that using a low (even though realistic) SNR might wash out the '
                           'higher frequencies, which could make the map unsuitable to be used for further refinement.')
        self._defineExtraParams(form)

    # -------------------------- INSERT steps functions -----------------------
    def _insertAllSteps(self):
        self._insertFunctionStep(self.convertInputStep, needsGPU=False)
        self._insertFunctionStep(self.relionReconstructParticle, needsGPU=False)
        self._insertFunctionStep(self.createOutputStep, needsGPU=False)

    # -------------------------- STEPS functions ------------------------------
    def relionReconstructParticle(self):
        cmd = self._genRecParticleCmd()
        try:
            Plugin.runRelionTomo(self, 'relion_tomo_reconstruct_particle_mpi', cmd,
                                 numberOfMpi=self.numberOfMpi.get())
        except:
            # The --mem argument should also be set using around 80-90% to keep a safety margin
            Plugin.runRelionTomo(self, 'relion_tomo_reconstruct_particle_mpi', cmd + '--mem 50 ',
                                 numberOfMpi=self.numberOfMpi.get())

    def createOutputStep(self):
        inParticles = self.inReParticles.get()
        currentSamplingRate = inParticles.getTsSamplingRate() * self.binningFactor.get()
        halves = [self._getExtraPath('half1.mrc'), self._getExtraPath('half2.mrc')]

        # Fix headers to be interpreted as volumes instead of stacks
        [fixVolume(mrcFile) for mrcFile in glob.glob(self._getExtraPath('*.mrc'))]

        # Output average
        vol = AverageSubTomogram()
        vol.setFileName(self._getExtraPath('merged.mrc'))
        if exists(halves[0]):
            vol.setHalfMaps(halves)
        vol.setSamplingRate(currentSamplingRate)
        self._defineOutputs(**{outputObjects.average.name: vol})
        self._defineSourceRelation(inParticles, vol)

    # -------------------------- INFO functions -------------------------------
    def _validate(self):
        validateMsg = super()._validate()
        if not validateMsg:
            if self.numberOfMpi.get() == 1:
                validateMsg.append('The number of MPI must be greater than 1')
        return validateMsg

    # --------------------------- UTILS functions -----------------------------
    def _genRecParticleCmd(self):
        cmd = [
            self._genCommonExtractAndRecCmd(),
            f'--sym {self.symmetry.get()}',
            # Note:
            #   --j: number of threads used for the non-reconstruction parts of the program (e.g. symmetry application
            #        or gridding correction). This should be set to the number of CPU cores available.
            #   --j_out: number of threads that compute partial reconstructions in parallel. This is faster, but it
            #        requires additional memory for each thread. When used together with the --mem argument, this number
            #        will be reduced to (approximately) maintain the imposed memory limitation.
            #   --j_in: number of threads to be used for each partial reconstruction. This is a slower way to
            #        parallelise the procedure, but it does not require additional memory. Unless memory is limited,
            #        the --j_out option should be preferred. The product of --j_out and --j_in should not exceed the
            #        number of CPU cores available.
            f'--j_out {self.binThreads.get()}',
            '--j_in 1'
        ]
        if self.snrWiener.get() > 0:
            cmd.append(f'--SNR {self.snrWiener.get():.2f}')
        return ' '.join(cmd)

