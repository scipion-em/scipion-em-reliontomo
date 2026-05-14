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
from os.path import exists, join
from typing import Union

from emtable import Table

from pwem.objects import Volume, FSC
from pwem.protocols import EMProtocol
from pyworkflow import BETA, PROD
from pyworkflow.object import Pointer
from pyworkflow.protocol import PointerParam, StringParam, IntParam
from pyworkflow.utils import Message, createLink
from reliontomo import Plugin
from reliontomo.constants import IN_PARTICLES_STAR, POSTPROCESS_DIR, OPTIMISATION_SET_STAR, PSUBTOMOS_SQLITE, \
    OUT_PARTICLES_STAR, OUT_TOMOS_STAR, TRAJECTORIES_STAR
from reliontomo.convert import writeSetOfPseudoSubtomograms, readSetOfPseudoSubtomograms, convert50_tomo
from reliontomo.objects import RelionSetOfPseudoSubtomograms
from tomo.objects import SetOfCoordinates3D

IS_RELION_50 = Plugin.isRe50()


class ProtRelionTomoBase(EMProtocol):
    """
    Provides the foundational infrastructure for cryo-electron tomography workflows based on RELION
    pseudo-subtomograms. The protocol establishes the common mechanisms required for handling particle
    metadata, generating compatible STAR files, managing tomographic reconstruction outputs, and
    maintaining interoperability between different RELION tomography processing stages.

    AI Generated:

    Relion Tomography Base Protocol (ProtRelionTomoBase) — User Manual
        Overview

        The Relion Tomography Base Protocol serves as the common foundation for multiple cryo-electron
        tomography workflows that operate on pseudo-subtomograms within the RELION environment. Its
        purpose is to standardize how tomographic particle information, metadata, reconstruction
        products, and intermediate processing files are handled throughout a tomography processing
        pipeline.

        In RELION tomography workflows, pseudo-subtomograms are used as an efficient representation of
        tomographic particle information. Rather than representing a direct physical reconstruction of
        the particle scattering potential, these objects provide a computationally practical framework
        for integrating two-dimensional tilt-series information into iterative refinement and averaging
        procedures. This approach allows RELION tomography to leverage established refinement strategies
        while remaining computationally tractable for large cryo-EM datasets.

        Pseudo-Subtomograms and Their Biological Role

        Pseudo-subtomograms are central to modern RELION tomography workflows because they provide the
        information required for alignment, refinement, and classification while preserving the
        relationship between particles and their original tilt-series observations. These objects are
        generated from contrast-transfer-function-corrected projections together with auxiliary
        information describing sampling statistics and frequency-space coverage.

        From a biological perspective, pseudo-subtomograms enable iterative refinement of macromolecular
        complexes observed inside cells, vesicles, viral assemblies, or purified specimens embedded in
        vitreous ice. They support workflows aimed at recovering higher-resolution structures while
        accounting for the geometric complexity inherent to tomography experiments.

        Input Data and Workflow Integration

        The protocol framework is designed to work with RELION pseudo-subtomogram datasets generated
        during earlier tomography processing stages. These datasets typically contain particle
        coordinates, alignment information, tomogram references, and associated metadata required for
        downstream refinement or reconstruction protocols.

        The protocol ensures that metadata is consistently transferred between processing stages,
        allowing subsequent workflows such as subtomogram averaging, classification, CTF refinement,
        polishing, or post-processing to access a coherent description of the experimental data.

        Compatibility Between RELION Versions

        An important role of the protocol is maintaining compatibility between tomography datasets
        generated using different RELION versions. RELION 4 and RELION 5 introduced differences in
        pseudo-subtomogram organization and metadata conventions, and this protocol framework provides
        mechanisms to validate dataset compatibility and generate the appropriate metadata structures
        required for each version.

        In practical terms, this helps users avoid inconsistencies that could otherwise lead to failed
        refinements, incorrect metadata interpretation, or incompatible processing pipelines. Users
        working across multiple computational environments or collaborating between facilities may
        particularly benefit from these compatibility safeguards.

        STAR File Management

        RELION tomography workflows rely heavily on STAR files to organize metadata describing particles,
        tomograms, trajectories, optics information, and refinement results. This protocol framework
        manages the generation and propagation of these files so that downstream protocols receive the
        required metadata in a consistent format.

        Proper STAR file management is biologically important because alignment parameters, particle
        identities, optical information, and reconstruction geometry must remain synchronized throughout
        the processing workflow. Even small metadata inconsistencies can compromise averaging quality or
        lead to incorrect structural interpretation.

        Reconstruction and Output Handling

        The protocol also provides mechanisms for handling tomographic reconstruction products,
        post-processed volumes, Fourier shell correlation measurements, and related metadata objects.
        These outputs are prepared in a form that allows direct integration into visualization,
        refinement, validation, and interpretation workflows within the Scipion ecosystem.

        Sampling rates, binning information, particle dimensions, and reconstruction geometry are
        propagated together with the reconstructed data to preserve physical consistency throughout the
        workflow. This is particularly important when combining multiple refinement stages or comparing
        reconstructions generated under different processing conditions.

        Computational Considerations

        Tomography processing can require substantial computational resources because many tilt-series
        projections, tomograms, and pseudo-subtomograms may need to be processed simultaneously. The
        protocol framework therefore supports configurable parallel execution strategies and integration
        with RELION multithreading capabilities.

        Users should consider available memory resources carefully when selecting processing parameters.
        Large tomograms, extensive particle sets, or highly parallel execution strategies can
        significantly increase memory usage and computational demand.

        Practical Recommendations

        For most biological workflows, maintaining metadata consistency across all tomography processing
        stages is essential. Users should ensure that pseudo-subtomogram datasets are generated using
        compatible RELION versions and that particle metadata remains synchronized after filtering,
        classification, or subset selection operations.

        When transferring datasets between workflows, it is advisable to verify sampling rates, box
        sizes, and tomography geometry before initiating downstream refinement or averaging steps.
        Maintaining consistent metadata throughout the pipeline reduces the risk of subtle alignment
        errors and improves the reproducibility of structural results.

        Final Perspective

        Modern cryo-electron tomography relies not only on accurate image processing algorithms but also
        on reliable management of complex metadata describing particles, optics, geometry, and
        reconstruction state. This protocol provides the organizational backbone required to maintain
        consistency throughout RELION tomography workflows, enabling robust subtomogram analysis and
        biologically meaningful structural interpretation.
    """
    _devStatus = BETA if IS_RELION_50 else PROD

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    @staticmethod
    def _defineCommonInputParams(form):
        form.addSection(label=Message.LABEL_INPUT)
        form.addParam('inReParticles', PointerParam,
                      important=True,
                      strict=True,
                      pointerClass='RelionSetOfPseudoSubtomograms',
                      label='Pseudo-Subtomograms',
                      help='Pseudo-subtomograms do not aim to accurately represent the scattering potential of '
                           'the underlying particles. Instead, they serve as a practical means to implement an '
                           'approximation to the 2D approach within the existing RELION framework. In the original '
                           'RELION4 article an accurate defition is given, see:\n '
                           'https://doi.org/10.7554/eLife.83724\n '
                           'A more technical explanation, pseudo-subtomograms are 3D-Arrays (volumes) that '
                           'constructed from the sums of 2D tilt-series images pre-multiplied by contrast transfer '
                           'functions (CTFs), along with auxiliary arrays that store the corresponding sum of squared '
                           'CTFs and the frequency of observation for each 3D voxel.')

    @staticmethod
    def _insertBinThreadsParam(form):
        form.addParam('binThreads', IntParam,
                      label='Relion threads',
                      important=True,
                      default=3,
                      help='Number of threads used by Relion each time it is called in the protocol execution. For '
                           'example, if 2 Scipion threads and 3 Relion threads are set, the tilt-series, tomograms, etc '
                           'will be processed in groups of 2 at the same time with a call of tomo3d with 3 threads '
                           'each, so 6 threads will be used at the same time. Beware the memory of your machine has '
                           'memory enough to load together the number of tilt-series, tomograms, etc specified by '
                           'Scipion threads.')

    @staticmethod
    def _defineExtraParams(form, addAdditionalSection=True):
        if addAdditionalSection:
            form.addSection(label='Additional')
        form.addParam('extraParams', StringParam,
                      label='Additional arguments',
                      help="In this box command-line arguments may be provided that are not generated by the GUI. This "
                           "may be useful for testing developmental options and/or expert use of the program, e.g: \n"
                           "--verb 1\n"
                           "--pad 2\n")

    # --------------------------- UTILS functions -----------------------------
    def getInputParticles(self, returnPointer: bool=False) -> Union[Pointer, RelionSetOfPseudoSubtomograms]:
        reParticlesPointer = self.inReParticles
        return reParticlesPointer if returnPointer else reParticlesPointer.get()

    def getOutStarFileName(self):
        return self._getExtraPath(IN_PARTICLES_STAR)

    def genInStarFile(self, withPriors=False, are2dParticles=False):
        """It will check if the set size and the stored particles star file are of the same size or not. In
        the first case, a link will be made to the previous particles star file to avoid generating it and in the
        second case, a new file will be generated containing only the ones present in the input set.
        :param withPriors: Applies only if using Relion 4. Consider the prior angles.
        :param are2dParticles: Applies only if using Relion 5. Used to choose the fields that will be present in
        the generated particles.star file, as they are not the same depending on if the particles are 2D or 3D.
        """
        inReParticlesSet = self.getInputParticles()
        outStarFileName = self.getOutStarFileName()
        # if IS_RELION_50:
        #     withPriors = False
        # if inReParticlesSet.getSize() == inReParticlesSet.getNReParticles() and not withPriors:
        #     if inReParticlesSet.getSize() == inReParticlesSet.getNReParticles():
        #         self.info("Using existing star (%s) file instead of generating a new one." %
        #                   inReParticlesSet.getParticlesStar())
        #         createLink(inReParticlesSet.getParticlesStar(), outStarFileName)
        # else:
        #     self.info("Less particles detected in the input set respecting to it associated star file. Assuming "
        #               "that a subset was made. Writing the new particles file.")
        if IS_RELION_50:
            outPath = self._getExtraPath()
            writer = convert50_tomo.Writer()
            writer.pseudoSubtomograms2Star(inReParticlesSet, outPath, are2dParticles=are2dParticles)
        else:
            writeSetOfPseudoSubtomograms(inReParticlesSet, outStarFileName, withPriors=withPriors)

    def _genPostProcessOutputMrcFile(self, fileName):
        """File generated using the sharpening protocol (called post-process protocol) and also using the
        rec particle from TS protocol in case the optional input 'solvent mask' is introduced."""
        postProccesMrc = Volume()
        postProccesMrc.setFileName(self._getExtraPath(POSTPROCESS_DIR, fileName))
        sRate = -1
        if getattr(self, 'inReParticles', None):
            sRate = self.inReParticles.get().getCurrentSamplingRate()
        elif getattr(self, 'inVolume', None):
            sRate = self.inVolume.get().getSamplingRate()
        postProccesMrc.setSamplingRate(sRate)

        return postProccesMrc

    def genRelionParticles(self,
                           optimisationFileName=OPTIMISATION_SET_STAR,
                           particles=OUT_PARTICLES_STAR,
                           tomograms=OUT_TOMOS_STAR,
                           trajectories=TRAJECTORIES_STAR,
                           binningFactor=None,
                           boxSize=24):
        """Generate a RelionSetOfPseudoSubtomograms object containing the files involved for the next protocol,
        considering that some protocols don't generate the optimisation_set.star file. In that case, the input Object
        which represents it will be copied and, after that, this method will be used to update the corresponding
        attribute."""

        # Create the set
        inParticlesSet = self.getInputParticles()
        psubtomoSet = RelionSetOfPseudoSubtomograms.create(self.getPath(), template=PSUBTOMOS_SQLITE)
        psubtomoSet.copyInfo(inParticlesSet)

        # Verify out star file
        optimSetStar =  self._getExtraPath(optimisationFileName)
        if exists(optimSetStar):
            psubtomoSet.filesMaster = optimSetStar

        particles =  self._getExtraPath(particles)
        if exists(particles):
            psubtomoSet.setParticles(particles)

        tomograms =  self._getExtraPath(tomograms)
        if exists(tomograms):
            psubtomoSet.setTomogramsStar(tomograms)

        trajectiories = self._getExtraPath(trajectories)
        if exists(trajectiories):
            psubtomoSet.setTrajectoriesStar(trajectories)

        if binningFactor:
            psubtomoSet.setRelionBinning(binningFactor)

        if boxSize:
            psubtomoSet.setBoxSize(boxSize)

        # Fill the items (pseudo subtomos/particles) from de particles star file
        psubtomoSet.setSamplingRate(psubtomoSet.getCurrentSamplingRate())
        readSetOfPseudoSubtomograms(psubtomoSet, isRelion5=IS_RELION_50)

        return psubtomoSet

    def genFSCs(self, starFile, tableName, fscColumns):
        fscSet = self._createSetOfFSCs()
        table = Table(fileName=starFile, tableName=tableName)
        resolution_inv = table.getColumnValues('rlnResolution')
        for columnName in fscColumns:
            columValues = table.getColumnValues(columnName)
            fsc = FSC(objLabel=columnName[3:])
            fsc.setData(resolution_inv, columValues)
            fscSet.append(fsc)

        fscSet.write()
        return fscSet

    def _genExtraParamsCmd(self):
        return ' ' + self.extraParams.get() if self.extraParams.get() else ''

    def isInputSetOf3dCoords(self):
        return True if type(self.getInputParticles()) is SetOfCoordinates3D else False

    # -------------------------- INFO functions -------------------------------
    def _validate(self):
        errorMsg = []
        inParticles = self.getInputParticles()
        if not self.isInputSetOf3dCoords():
            areRe5Particles = inParticles.areRe5Particles()
            if IS_RELION_50 and not areRe5Particles:
                errorMsg.append('The introduced particles were not generated with Relion 5, while the plugin is '
                                'currently configured to work with Relion 5. Please consider:'
                                '\n - Calling the protocol "Extract subtomos" to convert the particles into Relion 5 '
                                'format OR'
                                '\n - Configuring the plugin to work with Relion 4.')
            if not IS_RELION_50 and areRe5Particles:
                errorMsg.append('The introduced particles were generated with Relion 5, while the plugin is currently '
                                'configured to work with Relion 4. Please consider configuring the plugin to work with '
                                'Relion 5.')
        return errorMsg
