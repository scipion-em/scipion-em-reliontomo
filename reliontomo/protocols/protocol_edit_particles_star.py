# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors:     Scipion Team
# *
# * National Center of Biotechnology, CSIC, Spain
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
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************
from enum import Enum

from pyworkflow.object import Boolean
from pyworkflow.protocol import BooleanParam, FloatParam, EnumParam, \
    PointerParam, IntParam, GE
from reliontomo import Plugin
from reliontomo.constants import OUT_PARTICLES_STAR, COORD_X, COORD_Y, COORD_Z, SHIFTX_ANGST, SHIFTY_ANGST, \
    SHIFTZ_ANGST, ROT, TILT, PSI
from reliontomo.objects import RelionSetOfPseudoSubtomograms
from reliontomo.protocols.protocol_base_relion import ProtRelionTomoBase, IS_RELION_50
from reliontomo.utils import genEnumParamDict

# Operation labels and values
NO_OPERATION = 'No operation'
OP_ADDITION = 'Add'
OP_MULTIPLICATION = 'Multiply'
OP_SET_TO = 'Set to'
OPERATION_LABELS = [NO_OPERATION, OP_ADDITION, OP_MULTIPLICATION, OP_SET_TO]

# Labels to which apply the operation
COORDINATES = 'coordinates'
SHIFTS = 'shifts'
ANGLES = 'angles'
LABELS_TO_OPERATE_WITH = [COORDINATES, SHIFTS, ANGLES]


class outputObjects(Enum):
    relionParticles = RelionSetOfPseudoSubtomograms


class ProtRelionEditParticlesStar(ProtRelionTomoBase):
    """
    Performs editing, refinement, and cleanup operations on RELION particle metadata associated with subtomograms
    and tomography workflows. The protocol is intended to help users reorganize particle positions and orientations
    after alignment, classification, or refinement procedures in order to improve the consistency and biological
    interpretability of downstream analyses.

    AI Generated:

    Edit RELION Particles (ProtRelionEditParticlesStar) - User Manual
        Overview

        The Edit RELION Particles protocol provides a collection of utilities for modifying particle metadata
        associated with subtomogram averaging and cryo-electron tomography workflows. Its primary purpose is to
        adjust particle positioning, redefine particle centers, remove duplicated particles, and optionally modify
        geometric metadata such as coordinates, shifts, or angular assignments.

        In practical cryo-ET workflows, particle metadata often evolves throughout iterative alignment and refinement
        procedures. During these stages, particles may drift from their original centers, converge toward overlapping
        positions, or accumulate orientation adjustments that require correction before further refinement. This
        protocol serves as an intermediate curation and cleanup step that improves the reliability of subsequent
        analyses.

        Biological Context and Motivation

        In subtomogram averaging, accurate particle positioning is critical because even small centering errors can
        significantly reduce achievable resolution. Mis-centered particles introduce blurring during averaging and may
        obscure biologically meaningful structural details. For this reason, re-centering particles around a specific
        structural feature is often an essential refinement strategy.

        Duplicate particles represent another important issue in tomography workflows. During iterative refinement,
        neighboring particles may converge toward the same physical location, particularly in crowded cellular
        environments or membrane-associated assemblies. If not corrected, duplicated particles artificially inflate
        particle counts and compromise independent half-set validation procedures such as gold-standard FSC
        calculations.

        The protocol therefore supports operations that improve both the geometric consistency of the dataset and the
        statistical validity of downstream refinements.

        Inputs and General Workflow

        The protocol operates on RELION particle metadata generated from subtomogram averaging workflows. These
        particles typically contain positional information, angular assignments, translational shifts, and associated
        tomographic metadata required for downstream RELION processing.

        Users may optionally provide an average subtomogram and a reference mask when performing re-centering
        operations. These references help define the biologically relevant feature that should become the new particle
        center. This is particularly useful when the dominant signal in the particle is offset from the original
        extraction center or when refinement converges toward a distinct structural region.

        After processing, the protocol generates a new edited particle set that preserves the original particle
        identities while updating their geometric metadata according to the selected operations.

        Re-centering Particles

        Re-centering is one of the most biologically meaningful operations provided by this protocol. In many cryo-ET
        studies, the initial extraction coordinates correspond only approximately to the true structural center of the
        the complex. Subsequent refinement may reveal that a different region of the structure provides a more stable
        or biologically relevant alignment target.

        By shifting the particle center toward a chosen reference position, users can improve alignment consistency
        and reduce variability across the dataset. This is especially valuable for asymmetric complexes, elongated
        assemblies, membrane proteins, or particles containing flexible peripheral domains.

        The provided shifts correspond to positional offsets within the particle reference frame. Care should be taken
        to ensure that the chosen center reflects a structurally stable region rather than a flexible or poorly
        resolved component.

        Duplicate Particle Removal

        The protocol also supports removal of particles that occupy nearly identical spatial positions. This operation
        is particularly important in high-density cellular datasets where repeated refinement cycles may cause nearby
        particles to collapse toward the same coordinate.

        Removing duplicates improves dataset integrity and prevents overrepresentation of certain structural regions.
        It also protects the validity of resolution estimation procedures that rely on statistical independence
        between particle subsets.

        The minimum inter-particle distance parameter determines how aggressively duplicates are removed. Smaller
        values preserve closely packed particles, whereas larger values favor stricter cleanup. Biologically crowded
        systems such as ribosome-rich cytoplasm or membrane lattices may require conservative thresholds to avoid
        eliminating genuinely distinct particles.

        Metadata Editing Operations

        In earlier workflow configurations, the protocol additionally supports arithmetic operations on particle
        metadata fields such as coordinates, shifts, and Euler angles. These operations allow users to uniformly
        modify particle geometry across an entire dataset.

        Such functionality can be useful when correcting systematic offsets, compensating for known coordinate
        transformations, or adjusting orientation conventions between processing pipelines. However, these operations
        should be used cautiously because inappropriate modifications may disrupt the physical consistency of the
        dataset.

        Coordinate operations affect particle positions within the tomogram, shift operations modify alignment
        translations, and angular operations alter particle orientations. Since these parameters directly determine
        the spatial interpretation of the reconstruction, all modifications should be biologically justified and
        carefully validated visually.

        Outputs and Interpretation

        The protocol produces an updated RELION particle set suitable for subsequent subtomogram refinement,
        classification, averaging, or visualization workflows. The edited dataset maintains compatibility with RELION
        tomography processing while incorporating the selected geometric corrections.

        Re-centered particles typically yield improved alignment stability and sharper averages, whereas duplicate
        removal reduces statistical bias and improves reliability of resolution assessment. Metadata editing
        operations may also facilitate interoperability between different processing pipelines or experimental setups.

        Practical Recommendations

        In routine cryo-ET workflows, re-centering is most effective when guided by a well-resolved structural core.
        Users should avoid centering on highly flexible domains or poorly resolved peripheral regions, as this may
        destabilize refinement.

        Duplicate removal should generally be performed after major alignment refinements and before final resolution
        estimation. Conservative distance thresholds are recommended for crowded biological environments to avoid
        removing valid neighboring particles.

        Metadata editing operations should only be applied when there is a clear understanding of the coordinate and
        orientation conventions used throughout the workflow. Visual inspection of particle placement before and after
        editing is strongly recommended.

        Final Perspective

        Accurate particle geometry is fundamental for successful subtomogram averaging and structural interpretation.
        The Edit RELION Particles protocol provides a controlled environment for correcting positional inconsistencies,
        refining particle centering, and maintaining dataset quality throughout iterative tomography workflows.
        Careful use of these operations can substantially improve both reconstruction quality and biological
        interpretability.
    """

    _label = 'Apply operation to Relion particles'
    _possibleOutputs = outputObjects
    operationDict = genEnumParamDict(OPERATION_LABELS)
    labelsDict = genEnumParamDict(LABELS_TO_OPERATE_WITH)

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        super()._defineCommonInputParams(form)
        form.addSection(label='Center')
        form.addParam('doRecenter', BooleanParam,
                      label='Perform centering of particles',
                      default=False,
                      help='Perform centering of particles according to a position in the reference.')
        group = form.addGroup('Shift center', condition='doRecenter')
        group.addParam('averageSubTomogram',
                       PointerParam,
                       pointerClass='AverageSubTomogram',
                       label='Average of subtomogram (optional)',
                       allowsNull=True)
        group.addParam('refMask',
                       PointerParam,
                       pointerClass='VolumeMask',
                       label='Reference mask (optional)',
                       allowsNull=True)
        group.addParam('shiftX', FloatParam,
                       label='X (pix.)',
                       default=0,
                       help='X-coordinate in the reference to center particles on (in pix)')
        group.addParam('shiftY', FloatParam,
                       label='Y (pix.)',
                       default=0,
                       help='Y-coordinate in the reference to center particles on (in pix)')
        group.addParam('shiftZ', FloatParam,
                       label='Z (pix.)',
                       default=0,
                       help='Z-coordinate in the reference to center particles on (in pix)')
        # Section 'Remove duplicates' was added in Relion 5
        if IS_RELION_50:
            form.addSection('Remove duplicates')
            form.addParam('doRemoveDuplicates', BooleanParam,
                          default=False,
                          label='Remove duplicates?',
                          help='If set to Yes, duplicated particles that are within a given distance are removed '
                               'leaving only one. Duplicated particles are sometimes generated when particles drift '
                               'into the same position during alignment. They inflate and invalidate gold-standard '
                               'FSC calculation.')
            form.addParam('minDistPartRemoval', IntParam,
                          default=30,
                          label='Minimum inter-particle distance (Å)',
                          condition='doRemoveDuplicates',
                          validators=[GE(0)],
                          help='Particles within this distance are removed leaving only one.')
        else:
            # This section is rarely used and makes nonsense in Relion 5 as there is not a explicit field for the
            # coordinates
            form.addSection(label='Operate')
            form.addParam('chosenOperation', EnumParam,
                          choices=list(self.operationDict.keys()),
                          default=self.operationDict[NO_OPERATION],
                          label='Choose operation')
            form.addParam('opValue', FloatParam,
                          condition='chosenOperation != %i' % self.operationDict[NO_OPERATION],
                          default=1,
                          label='Value to operate the selected labels')
            group = form.addGroup('Operation', condition='chosenOperation > 0')
            group.addParam('operateWith', EnumParam,
                           choices=list(self.labelsDict.keys()),
                           default=self.labelsDict[COORDINATES],
                           label='Operate with')
            group.addParam('label1x', BooleanParam,
                           label='X (pix.)',
                           condition='operateWith == %s' % self.labelsDict[COORDINATES],
                           default=False)
            group.addParam('label2y', BooleanParam,
                           label='Y (pix.)',
                           condition='operateWith == %s' % self.labelsDict[COORDINATES],
                           default=False)
            group.addParam('label3z', BooleanParam,
                           label='Z (pix.)',
                           condition='operateWith == %s' % self.labelsDict[COORDINATES],
                           default=False)
            group.addParam('label1sx', BooleanParam,
                           label='Shift X (pix.)',
                           condition='operateWith == %s' % self.labelsDict[SHIFTS],
                           default=False)
            group.addParam('label2sy', BooleanParam,
                           label='Shift Y (pix.)',
                           condition='operateWith == %s' % self.labelsDict[SHIFTS],
                           default=False)
            group.addParam('label3sz', BooleanParam,
                           label='Shift Z (pix.)',
                           condition='operateWith == %s' % self.labelsDict[SHIFTS],
                           default=False)
            group.addParam('label1rot', BooleanParam,
                           label='Rot (deg.)',
                           condition='operateWith == %s' % self.labelsDict[ANGLES],
                           default=False)
            group.addParam('label2tilt', BooleanParam,
                           label='Tilt (deg.)',
                           condition='operateWith == %s' % self.labelsDict[ANGLES],
                           default=False)
            group.addParam('label3psi', BooleanParam,
                           label='Psi (deg.)',
                           condition='operateWith == %s' % self.labelsDict[ANGLES],
                           default=False)

    # -------------------------- INSERT steps functions -----------------------
    def _insertAllSteps(self):
        self._insertFunctionStep(self.convertInputStep, needsGPU=False)
        self._insertFunctionStep(self.operateStep, needsGPU=False)
        self._insertFunctionStep(self.createOutputStep, needsGPU=False)

    def convertInputStep(self):
        are2dParticles = self.inReParticles.get().areRe5Particles()
        self.genInStarFile(are2dParticles=are2dParticles)

    def operateStep(self):
        program = 'relion_star_handler'
        if IS_RELION_50:
            # Remove duplicates and re-center particles needs to be executed separately, being the second call to the
            # program required to be fed with the ourpur star of the first
            doRemDuplicates = self.doRemoveDuplicates.get()
            doRecenter = self.doRecenter.get()
            if doRemDuplicates and doRecenter:
                Plugin.runRelionTomo(self, program, self._genRemDuplicatesCmd())
                outParticles = self._getExtraPath(OUT_PARTICLES_STAR)
                cmd = '--i %s ' % outParticles
                cmd += '--o %s ' % outParticles
                cmd += self._genRecenterCmd()
                Plugin.runRelionTomo(self, program, cmd)
            else:
                if doRemDuplicates:
                    Plugin.runRelionTomo(self, program, self._genRemDuplicatesCmd())
                if doRecenter:
                    Plugin.runRelionTomo(self, program, self._genRecenterCmdWithIO())
        else:
            Plugin.runRelionTomo(self, program, self._getOperateCommand())

    def createOutputStep(self):
        psubtomoSet = self.genRelionParticles()
        self._defineOutputs(**{outputObjects.relionParticles.name: psubtomoSet})
        self._defineSourceRelation(self.inReParticles, psubtomoSet)

    # -------------------------- INFO functions -------------------------------
    def _validate(self):
        valMsg = []
        if IS_RELION_50:
            if not self.doRecenter.get() and not self.doRemoveDuplicates.get():
                valMsg.append('No re-centering nor duplicate removal operation was chosen.')
        else:
            if not self.doRecenter.get() and self.chosenOperation.get() == self.operationDict[NO_OPERATION]:
                valMsg.append('No re-centering or operation was chosen.')
        return valMsg

    # --------------------------- UTILS functions -----------------------------
    def _genIOCmd(self):
        cmd = '--i %s ' % self.getOutStarFileName()
        cmd += '--o %s ' % self._getExtraPath(OUT_PARTICLES_STAR)
        return cmd

    def _genRemDuplicatesCmd(self):
        cmd = self._genIOCmd()
        cmd += '--remove_duplicates %i ' % self.minDistPartRemoval.get()
        return cmd

    def _genRecenterCmd(self):
        cmd = '--center '
        if self.shiftX.get() != 0:
            cmd += '--center_X %.2f ' % self.shiftX.get()
        if self.shiftY.get() != 0:
            cmd += '--center_Y %.2f ' % self.shiftY.get()
        if self.shiftZ.get() != 0:
            cmd += '--center_Z %.2f ' % self.shiftZ.get()
        return cmd

    def _genRecenterCmdWithIO(self):
        cmd = self._genIOCmd()
        cmd += self._genRecenterCmd()
        return cmd

    def _getOperateCommand(self):
        cmd = self._genIOCmd()
        # Re-center particles
        if self.doRecenter.get():
            cmd += self._genRecenterCmd()

        # Operate particles - removed by Jorge in the protocol version for Relion 5 as it is rarely used and may
        # be problematic as some of the fields involved are now defined in a different way, such as the coordinates
        if self.chosenOperation.get() != self.operationDict[NO_OPERATION]:
            opValue = self.opValue.get()
            chosenOp = self.chosenOperation.get()
            # Chosen operation
            if chosenOp == self.operationDict[OP_ADDITION]:
                cmd += '--add_to %.2f ' % opValue
            elif chosenOp == self.operationDict[OP_MULTIPLICATION]:
                cmd += '--multiply_by %.2f ' % opValue
            else:
                cmd += '--set_to %.2f ' % opValue
            # Chosen values
            cmd += self._genOperateCmd()
        return cmd

    def _genOperateCmd(self):
        """Three are the maximum number of labels able to be edited at once. Relion offers 3 arguments to add them
        to the generated command: --operate, --operate2, --operate3."""
        operateWith = self.operateWith.get()
        if operateWith == self.labelsDict[COORDINATES]:
            label1, label2, label3 = COORD_X, COORD_Y, COORD_Z
            edit1, edit2, edit3 = self.label1x.get(), self.label2y.get(), self.label3z.get()
        elif operateWith == self.labelsDict[ANGLES]:
            label1, label2, label3 = ROT, TILT, PSI
            edit1, edit2, edit3 = self.label1rot.get(), self.label2tilt.get(), self.label3psi.get()
        else:
            label1, label2, label3 = SHIFTX_ANGST, SHIFTY_ANGST, SHIFTZ_ANGST
            edit1, edit2, edit3 = self.label1sx.get(), self.label2sy.get(), self.label3sz.get()

        operateCmd = ''
        labelList = [label1, label2, label3]
        editList = [edit1, edit2, edit3]
        counter = 1
        for label, editVal in zip(labelList, editList):
            if editVal:
                if counter == 1:
                    operateCmd += '--operate %s ' % label
                else:
                    operateCmd += '--operate%i %s ' % (counter, label)
                counter += 1

        return operateCmd

    def getAverageSubTomogram(self):
        return self.averageSubTomogram.get()

    def getMask3D(self):
        return self.refMask.get()
