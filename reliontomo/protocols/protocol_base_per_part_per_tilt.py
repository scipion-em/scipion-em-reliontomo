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
import shutil
from enum import Enum
from os import walk
from os.path import join, basename, exists
from typing import List
from emtable import Table
from pyworkflow.object import String
from pyworkflow.protocol import PointerParam, IntParam, GE
from pyworkflow.utils import createLink
from reliontomo.constants import IN_TOMOS_STAR, OUT_TOMOS_STAR, GLOBAL_TABLE, tomoStarFields, RLN_TOMONAME, RLN_VOLTAGE, \
    RLN_SPHERICALABERRATION, RLN_AMPLITUDECONTRAST, RLN_MICROGRAPHORIGINALPIXELSIZE, RLN_TOMOHAND, RLN_OPTICSGROUPNAME, \
    RLN_TOMOTILT_SERIES_PIXEL_SIZE, RLN_ETOMO_DIRECTIVE_FILE, RLN_TOMOTOMOGRAM_BINNING, RLN_TOMOSIZEX, RLN_TOMOSIZEY, \
    RLN_TOMOSIZEZ, RLN_TOMORECONSTRUCTED_TOMOGRAM, POSTPROCESS_STAR_FIELD
from reliontomo.objects import RelionSetOfPseudoSubtomograms
from reliontomo.protocols.protocol_base_relion import ProtRelionTomoBase, IS_RELION_50


class outputObjects(Enum):
    relionParticles = RelionSetOfPseudoSubtomograms


class ProtRelionPerParticlePerTiltBase(ProtRelionTomoBase):
    """
    Base protocol used for the getting the frame alignment and ctf-refinment

    AI Generated:

    Per-Particle and Per-Tilt Refinement Base
    (ProtRelionPerParticlePerTiltBase) - User Manual

        Overview

        This protocol provides the common framework for advanced
        tomography refinement procedures that operate at the level of
        individual particles and individual tilt images. Its main goal
        is to improve the quality of pseudo-subtomograms by refining
        local image properties that cannot be fully corrected during
        earlier reconstruction stages. These procedures are especially
        important in high-resolution cryo-electron tomography workflows,
        where small alignment or optical inaccuracies can strongly
        affect the final structural interpretation.

        In biological practice, this type of refinement is commonly
        used after an initial reconstruction has already been obtained.
        The protocol assumes that particles have already been extracted,
        aligned, and associated with tomographic metadata. It then uses
        a refined reference together with tomographic geometry
        information to improve consistency between the experimental
        tilt-series data and the reconstructed structure.

        Biological Purpose

        Cryo-electron tomography datasets often contain local motion,
        beam-induced distortions, or residual optical inaccuracies that
        vary from tilt to tilt and from particle to particle. Global
        corrections alone are frequently insufficient for demanding
        structural studies. This protocol establishes the infrastructure
        needed for more precise local refinement procedures that improve
        map sharpness, alignment stability, and high-resolution signal
        recovery.

        From a biological perspective, these refinements become
        especially important for heterogeneous specimens, flexible
        assemblies, membrane-associated complexes, and crowded cellular
        environments. Subtle structural details that are initially
        obscured by local inaccuracies may become significantly clearer
        after refinement.

        Inputs and Reference Information

        The protocol requires pseudo-subtomograms together with their
        associated tomographic metadata. These particles must already be
        connected to the original tilt-series information so that local
        corrections can be propagated consistently through the workflow.

        A reconstructed volume with associated half maps is also
        required. The half maps are essential because they provide an
        independent estimation of signal reliability and reduce the risk
        of overfitting during refinement. In practical workflows, the
        reference should ideally correspond to the best available
        reconstruction of the biological state under study.

        An optional soft reference mask can be provided to focus the
        refinement on biologically meaningful regions. This is
        particularly useful for flexible complexes, membrane proteins,
        or assemblies embedded within noisy tomographic environments.
        Proper masking helps emphasize stable structural regions while
        reducing the influence of surrounding solvent or unrelated
        density.

        FSC Information and Signal Reliability

        The protocol optionally accepts Fourier Shell Correlation
        information derived from post-processing procedures. This data
        improves the estimation of signal-to-noise relationships during
        refinement and provides a more realistic interpretation of
        resolution-dependent confidence.

        When FSC information is unavailable, refinement can still
        proceed, but the estimated signal quality may become slightly
        optimistic. For exploratory analyses this is often acceptable,
        whereas for high-resolution or publication-level studies the use
        of properly validated FSC information is strongly recommended.

        Box Size and High-Resolution Recovery

        The estimation box size defines the spatial region used during
        refinement calculations. Choosing an appropriate value is
        biologically important because high-frequency information may be
        spatially delocalized by the contrast transfer function.

        In practice, larger estimation boxes can improve the recovery
        of high-resolution information, particularly for large
        complexes, thick specimens, or datasets collected under strong
        defocus conditions. However, larger boxes also increase memory
        usage and computational cost. A balanced choice is therefore
        recommended depending on specimen size and desired resolution.

        Workflow Integration

        This protocol is designed as a foundational component for more
        specialized refinement procedures such as motion correction,
        trajectory refinement, or CTF refinement within tomographic
        workflows. It manages the relationships between particle
        metadata, tomographic geometry, reference information, and
        refinement outputs so that downstream protocols can operate
        consistently.

        During execution, the protocol preserves and updates the
        tomographic metadata associated with each particle. This ensures
        that subsequent processing stages continue to operate on a
        coherent representation of the experimental dataset.

        Outputs and Their Interpretation

        The main output is an updated set of pseudo-subtomograms with
        refined metadata and improved consistency relative to the
        reconstructed reference. These refined particles are intended
        for subsequent high-resolution refinement, classification, or
        reconstruction procedures.

        In practical biological workflows, improvements may appear as
        enhanced map sharpness, better-resolved secondary structure,
        clearer membrane features, or improved interpretability of
        flexible regions. The quality of these improvements strongly
        depends on the accuracy of the initial reconstruction and the
        consistency of the underlying tomographic data.

        Practical Recommendations

        For most datasets, the reference volume should already be well
        refined before applying per-particle or per-tilt corrections.
        Poor initial references may propagate incorrect structural
        features into later stages of refinement.

        The reference mask should include the biologically meaningful
        signal while excluding excessive solvent regions. Overly tight
        masks may suppress flexible but important structural features,
        whereas overly loose masks may reduce refinement sensitivity.

        Providing validated FSC information is recommended whenever
        reliable resolution estimation is important. This is especially
        relevant in high-resolution subtomogram averaging projects.

        Because these procedures can become computationally demanding,
        users should also ensure that sufficient processing resources
        and storage capacity are available before starting large-scale
        refinements.

        Final Perspective

        Per-particle and per-tilt refinement represents one of the most
        important advances in modern cryo-electron tomography because it
        enables correction of local inaccuracies that limit structural
        resolution. By combining particle metadata, tilt-series
        geometry, half-map validation, and localized refinement
        strategies, this protocol establishes the foundation for
        obtaining biologically meaningful high-resolution tomographic
        reconstructions.
    """

    _possibleOutputs = outputObjects

    # -------------------------- DEFINE param functions -----------------------

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.inParticlesStar = None

    def _defineParams(self, form):
        super()._defineCommonInputParams(form)
        self._insertBinThreadsParam(form)
        form.addParam('recVolume', PointerParam,
                      pointerClass='AverageSubTomogram',
                      allowsNull=False,
                      important=True,
                      label='Volume to get the halves',
                      help='Provide a volume with half maps. Note that volumes with associated'
                           'halves in the Scipion summary as w/h.')
        form.addParam('inRefMask', PointerParam,
                      pointerClass='VolumeMask',
                      label="Input reference mask",
                      help='This mask localizes the signal in the reference map. The mask should be'
                           'soft (non-binary)')
        if IS_RELION_50:
            form.addParam('inFsc', PointerParam,
                          pointerClass='SetOfFSCs',
                          allowsNull=True,
                          label='Post-process FSC (opt)',
                          help='If not provided, the SNR will be calculated without phase randomization, '
                               'so it will be slightly optimistic.')

    @staticmethod
    def _insertBoxSizeForEstimationParam(form):
        form.addParam('boxSize', IntParam,
                      label='Box size for estimation (px)',
                      default=128,
                      allowsNull=False,
                      validators=[GE(16)],
                      help="Box size to be used for the estimation. Note that this can be larger than the box size "
                           "of the reference map. A sufficiently large box size allows more of the high-frequency "
                           "signal to be captured that has been delocalized by the CTF.")

    # -------------------------- INSERT steps functions -----------------------
    def _insertAllSteps(self):
        # To be defined in the subclasses
        pass

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

    def createOutputStep(self):
        if IS_RELION_50:
            # Manage the output files
            tsStarDirTree = self._getExtraPath('Runs')
            resultingTsStarFiles = self._getResultingTsStarFiles(tsStarDirTree)
            outTsStarList = []
            for inTsStarFile in resultingTsStarFiles:
                outTsStar = self._getExtraPath(basename(inTsStarFile))
                shutil.move(inTsStarFile, outTsStar)
                outTsStarList.append(outTsStar)
            outTsStarList = sorted(outTsStarList)
            if exists(tsStarDirTree):
                shutil.rmtree(tsStarDirTree)
            # Update the file tomograms.star with the outTsStarFiles
            outTomoStar = self._getExtraPath(OUT_TOMOS_STAR)
            dataTable = Table()
            dataTable.read(outTomoStar, tableName=GLOBAL_TABLE)
            dataTable.sort(RLN_TOMONAME)
            tomoTable = Table(columns=tomoStarFields)
            with open(outTomoStar, 'w') as f:
                for row, tsStarFile in zip(dataTable, outTsStarList):
                    tomoTable.addRow(
                        row.get(RLN_TOMONAME),
                        row.get(RLN_VOLTAGE),
                        row.get(RLN_SPHERICALABERRATION),
                        row.get(RLN_AMPLITUDECONTRAST),
                        row.get(RLN_MICROGRAPHORIGINALPIXELSIZE),
                        row.get(RLN_TOMOHAND),
                        row.get(RLN_OPTICSGROUPNAME),
                        row.get(RLN_TOMOTILT_SERIES_PIXEL_SIZE),
                        tsStarFile,
                        row.get(RLN_ETOMO_DIRECTIVE_FILE),
                        row.get(RLN_TOMOTOMOGRAM_BINNING),
                        row.get(RLN_TOMOSIZEX),
                        row.get(RLN_TOMOSIZEY),
                        row.get(RLN_TOMOSIZEZ),
                        row.get(RLN_TOMORECONSTRUCTED_TOMOGRAM)
                    )
                # Write the STAR file
                tomoTable.writeStar(f, tableName=GLOBAL_TABLE)

        # Register outputs and define relations
        pSubtomos = self.genRelionParticles()  # Output RelionParticles
        self._defineOutputs(**{outputObjects.relionParticles.name: pSubtomos})
        self._defineSourceRelation(self.getInputParticles(), pSubtomos)

    # -------------------------- UTILS functions -----------------------------
    def _genIOCommand(self):
        inPSubtomos = self.inReParticles.get()
        inVolume = self.recVolume.get()
        trajectories = inPSubtomos.getTrajectoriesStar()
        postProcess = self._getPostProcessStar()
        half1, half2 = inVolume.getHalfMaps().split(',')
        cmd = '--p %s ' % inPSubtomos.getParticlesStar()
        cmd += '--t %s ' % inPSubtomos.getTomogramsStar()
        cmd += '--o %s ' % self._getExtraPath()
        if trajectories:
            cmd += '--mot %s ' % trajectories
        cmd += '--ref1 %s ' % half1
        cmd += '--ref2 %s ' % half2
        cmd += '--mask %s ' % self.inRefMask.get().getFileName()
        if postProcess:
            cmd += '--fsc %s ' % postProcess
        cmd += '--b %i ' % self.boxSize.get()
        cmd += '--j %i ' % self.binThreads.get()
        cmd += self._genExtraParamsCmd()
        return cmd

    @staticmethod
    def _getResultingTsStarFiles(iDir, ext='.star') -> List:
        tsStaFiles = []
        for dirpath, dirnames, filenames in walk(iDir):
            for filename in filenames:
                if filename.lower().endswith(ext):  # Ensure case-insensitive search
                    filepath = join(dirpath, filename)
                    tsStaFiles.append(filepath)
        return tsStaFiles

    def _getPostProcessStar(self):
        fsc = self.inFsc.get()
        if fsc:
            mdPostProcessStar = getattr(fsc, POSTPROCESS_STAR_FIELD, String()).get()
            if mdPostProcessStar and exists(mdPostProcessStar):
                return mdPostProcessStar
        return None

    # -------------------------- INFO functions -------------------------------
    def _warnings(self):
        warnMsg = []
        postProcessStar = self._getPostProcessStar()
        if not postProcessStar:
            warnMsg.append('No FSC was introduced or no postprocess.star file was found in the metadata of the '
                           'introduced FSC.\n'
                           'In that case, the SNR will be calculated without phase randomization, so it will be '
                           'slightly optimistic.')
        return warnMsg
