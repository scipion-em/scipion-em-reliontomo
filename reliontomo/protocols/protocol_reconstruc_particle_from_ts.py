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

from pwem.convert.headers import fixVolume
from pwem.objects import VolumeMask
from pyworkflow.protocol import StringParam, PointerParam, FloatParam
from reliontomo import Plugin
from reliontomo.constants import SYMMETRY_HELP_MSG, POST_PROCESS_MRC
from reliontomo.objects import RelionSetOfPseudoSubtomograms
from reliontomo.protocols.protocol_base_make_pseusosubtomos_and_rec_particle import \
    ProtRelionMakePseudoSubtomoAndRecParticleBase
from tomo.objects import AverageSubTomogram


class outputObjects(Enum):
    average = AverageSubTomogram
    postProcessVolume = VolumeMask
    relionParticles = RelionSetOfPseudoSubtomograms


class ProtRelionReconstructParticle(ProtRelionMakePseudoSubtomoAndRecParticleBase):
    """
    Reconstructs three-dimensional subtomogram averages from aligned tilt
    series particle projections. The protocol combines information from
    multiple tomographic views to generate refined volumetric
    reconstructions suitable for structural interpretation, downstream
    refinement, and subtomogram averaging workflows.

    AI Generated:

    Reconstruct Particle (ProtRelionReconstructParticle) —
    User Manual

        Overview

        This protocol reconstructs three-dimensional particle volumes from
        aligned cryo-electron tomography tilt-series data using Relion
        tomography reconstruction strategies. Its primary purpose is to
        generate averaged subtomogram maps that integrate information from
        multiple particle projections collected across the tilt series.

        In biological cryo-electron tomography workflows, this reconstruction
        step is essential because it transforms aligned particle information
        into interpretable three-dimensional density maps. These
        reconstructions often serve as the foundation for downstream
        refinement, classification, structural interpretation, or
        high-resolution subtomogram averaging.

        The protocol is particularly useful when generating consensus
        reconstructions from particles already aligned within tomographic
        datasets. By combining signal from many projections, the method
        improves the signal-to-noise ratio and enhances structural detail
        that may not be visible in individual tilt images.

        Inputs and Reconstruction Workflow

        The protocol requires aligned tomographic particle information derived
        from previous subtomogram alignment or pseudo-subtomogram generation
        workflows. These particles contain the positional and orientational
        information necessary to reconstruct a common three-dimensional map.

        During reconstruction, information from multiple tilt projections is
        merged into volumetric space to generate a consensus structural
        representation. The resulting map reflects the average structural
        organization of the particle population under analysis.

        In practical biological workflows, reconstruction quality depends
        strongly on the quality of the preceding alignment procedures.
        Accurate particle positioning and orientation assignments are critical
        because even small inconsistencies can blur high-resolution structural
        information in the final map.

        Symmetry Application

        The protocol allows the application of molecular symmetry during
        reconstruction. Correct symmetry specification can substantially
        improve reconstruction quality by averaging equivalent structural
        views and increasing the effective signal-to-noise ratio.

        Symmetry is particularly beneficial for highly ordered biological
        assemblies such as viral capsids, oligomeric protein complexes, or
        symmetric molecular machines. However, imposing incorrect symmetry may
        artificially distort asymmetric structural features or obscure
        biologically meaningful conformational variability.

        For flexible, asymmetric, or heterogeneous particles, reconstruction
        without symmetry constraints is generally the safest strategy. Users
        should select symmetry carefully according to the known biological
        organization of the specimen.

        Wiener Filtering and Signal Recovery

        The protocol supports Wiener filtering approaches that regulate the
        balance between signal preservation and noise suppression during
        reconstruction. This filtering strategy is biologically important
        because cryo-electron tomography data typically contain substantial
        levels of experimental noise.

        Appropriate filtering may improve visual interpretability and reduce
        reconstruction artifacts, especially in noisy or low-contrast
        datasets. However, excessively aggressive filtering can suppress
        biologically meaningful high-resolution information and may produce
        overly smooth reconstructions.

        In practice, moderate filtering values are often preferable when the
        reconstructed volume will be used for further refinement. Strong noise
        suppression may improve visualization but can reduce the ability of
        subsequent refinement procedures to recover high-frequency detail.

        Solvent Masking and FSC Estimation

        An optional solvent mask may be provided to support post-processing
        and Fourier Shell Correlation estimation. From a biological
        perspective, masking helps isolate the molecular region from the
        surrounding solvent background and improves the robustness of
        resolution assessment.

        Soft masks are generally preferred because they reduce edge artifacts
        and minimize artificial high-frequency effects introduced by sharp
        boundaries. Proper masking is especially important for irregularly
        shaped complexes, membrane proteins, or assemblies containing flexible
        peripheral regions.

        The protocol can also generate masked post-processed reconstructions
        suitable for downstream visualization and validation workflows.
        Resolution estimates derived from masked half-maps often provide more
        realistic assessments of structural quality than unmasked comparisons.

        Computational Considerations

        Three-dimensional subtomogram reconstruction is computationally
        demanding because large tomographic datasets require extensive memory
        and parallel processing resources. The protocol is designed to operate
        efficiently in distributed computing environments using MPI-based
        parallelization.

        Reconstruction performance depends on factors such as particle number,
        box size, symmetry, and reconstruction sampling. Larger particles and
        higher-resolution reconstructions generally require substantially more
        computational resources.

        In practical workflows, memory optimization strategies may become
        necessary when processing large datasets or high-resolution
        reconstructions. Computational efficiency should be balanced carefully
        against reconstruction quality and numerical stability.

        Outputs and Biological Interpretation

        The protocol produces reconstructed subtomogram averages together with
        independently reconstructed half-maps suitable for validation and
        resolution estimation. These outputs provide the basis for biological
        interpretation and downstream structural analysis.

        The averaged reconstruction represents the consensus structural state
        recovered from the particle population. Improved averaging often
        reveals secondary structure organization, domain architecture,
        membrane organization, or ligand-binding regions that are difficult
        to observe in raw tomographic data.

        Half-maps are particularly important because they allow independent
        validation of reconstruction quality through Fourier Shell Correlation
        analysis. These half-maps may also be used for subsequent
        post-processing procedures such as sharpening or local resolution
        estimation.

        Updated particle metadata may additionally be generated to support
        iterative refinement workflows. These refined particles can be reused
        in subsequent subtomogram averaging or classification procedures.

        Practical Recommendations

        In most biological workflows, reconstruction quality benefits greatly
        from accurate particle alignment and careful preprocessing of the
        tilt-series data. Users should verify alignment quality before
        attempting high-resolution reconstruction.

        Symmetry should only be applied when strongly supported by biological
        evidence. Incorrect symmetry assignment is one of the most common
        causes of reconstruction artifacts and misleading structural
        interpretation.

        Moderate filtering and conservative masking strategies are generally
        preferable during early reconstruction stages. More aggressive
        post-processing should typically be reserved for visualization and
        presentation purposes after structural validity has been confirmed.

        Final Perspective

        Subtomogram reconstruction is a central step in cryo-electron
        tomography because it transforms noisy projection data into coherent
        three-dimensional structural information. The quality of the final
        reconstruction strongly influences all downstream biological
        interpretation.

        Careful alignment, appropriate symmetry selection, thoughtful masking,
        and realistic noise treatment are essential for obtaining reliable and
        biologically meaningful reconstructions from tomographic datasets.
    """

    _label = 'Reconstruct particle'
    _possibleOutputs = outputObjects

    def __init__(self, **args):
        super().__init__(**args)

    # -------------------------- DEFINE param functions -----------------------

    def _defineParams(self, form):
        super()._defineParams(form)
        form.addSection(label='Reconstruct particle')
        super()._defineCommonRecParams(form)
        form.addParam('symmetry', StringParam,
                      label='Symmetry group',
                      default='C1',
                      help=SYMMETRY_HELP_MSG)
        form.addParam('solventMask', PointerParam,
                      pointerClass='VolumeMask',
                      label='FSC solvent mask (opt.)',
                      allowsNull=True,
                      help='Provide a soft mask to automatically estimate the postprocess FSC.')
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
        self._insertFunctionStep(self.convertInputStep)
        self._insertFunctionStep(self.relionReconstructParticle)
        if self.solventMask.get():
            self._insertFunctionStep(self.relionTomoMaskReference)
        self._insertFunctionStep(self.createOutputStep)

    # -------------------------- STEPS functions ------------------------------
    def relionReconstructParticle(self):
        cmd = self._genRecParticleCmd()
        try:
            Plugin.runRelionTomo(self, 'relion_tomo_reconstruct_particle_mpi', cmd, numberOfMpi=self.numberOfMpi.get())
        except:
            # The --mem argument should also be set using around 80-90% to keep a safety margin
            Plugin.runRelionTomo(self, 'relion_tomo_reconstruct_particle_mpi', cmd + '--mem 50 ',
                                 numberOfMpi=self.numberOfMpi.get())

    def relionTomoMaskReference(self):
        Plugin.runRelionTomo(self, 'relion_tomo_make_reference', self._genTomoMaskRefCmd(),
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
        vol.setHalfMaps(halves)
        vol.setSamplingRate(currentSamplingRate)
        self._defineOutputs(**{outputObjects.average.name: vol})

        # Output solvent mask
        if self.solventMask.get():
            postProccesMrc = self._genPostProcessOutputMrcFile(POST_PROCESS_MRC)
            postProccesMrc.setHalfMaps(halves)
            postProccesMrc.setSamplingRate(currentSamplingRate)
            self._defineOutputs(**{outputObjects.postProcessVolume.name: postProccesMrc})
            self._defineSourceRelation(inParticles, postProccesMrc)

        self._defineSourceRelation(inParticles, vol)

        # Create the output set with the new optimization set
        outParticles = self.genRelionParticles(boxSize=self.boxSize.get(),
                                               binningFactor=self.binningFactor.get())
        self._defineOutputs(**{outputObjects.relionParticles.name: outParticles})

    # -------------------------- INFO functions -------------------------------
    def _validate(self):
        validateMsg = []
        if self.numberOfMpi.get() == 1:
            validateMsg.append('The number of MPI must be greater than 1')
        return validateMsg

    # --------------------------- UTILS functions -----------------------------
    def _genRecParticleCmd(self):
        cmd = self._genCommonCmd()
        cmd += '--o %s ' % self._getExtraPath()
        cmd += '--sym %s ' % self.symmetry.get()
        # Note:
        #   --j: number of threads used for the non-reconstruction parts of the program (e.g. symmetry application
        #        or gridding correction). This should be set to the number of CPU cores available.
        #   --j_out: number of threads that compute partial reconstructions in parallel. This is faster, but it
        #        requires additional memory for each thread. When used together with the --mem argument, this number
        #        will be reduced to (approximately) maintain the imposed memory limitation.
        #   --j_in: number of threads to be used for each partial reconstruction. This is a slower way to parallelise
        #        the procedure, but it does not require additional memory. Unless memory is limited, the --j_out option
        #        should be preferred. The product of --j_out and --j_in should not exceed the number of CPU cores
        #        available.
        cmd += '--j_out %i ' % self.binThreads.get()
        cmd += '--j_in %i ' % 1
        if self.snrWiener.get() > 0:
            cmd += '--SNR %.2f ' % self.snrWiener.get()
        return cmd

    def _genTomoMaskRefCmd(self):
        inParticles = self.inReParticles.get()
        cmd = ''
        cmd += '--t %s ' % inParticles.getTomogramsStar()
        cmd += '--p %s ' % self.getOutStarFileName()
        cmd += '--rec %s ' % self._getExtraPath()
        cmd += '--o %s ' % self._getExtraPath()
        cmd += '--mask %s ' % self.solventMask.get().getFileName()
        cmd += '--angpix %.2f ' % (inParticles.getTsSamplingRate() * self.binningFactor.get())
        return cmd
