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
from pyworkflow.protocol import FloatParam, BooleanParam
from reliontomo import Plugin
from reliontomo.constants import OPTIMISATION_SET_STAR, PSUBTOMOS_SQLITE
from reliontomo.convert import writeSetOfPseudoSubtomograms, readSetOfPseudoSubtomograms
from reliontomo.objects import RelionSetOfPseudoSubtomograms, createSetOfRelionPSubtomograms
from reliontomo.protocols.protocol_base_make_pseusosubtomos_and_rec_particle import \
    ProtRelionMakePseudoSubtomoAndRecParticleBase
from reliontomo.utils import getProgram
from tomo.protocols import ProtTomoBase


class outputObjects(Enum):
    relionParticles = RelionSetOfPseudoSubtomograms


class ProtRelionMakePseudoSubtomograms(ProtRelionMakePseudoSubtomoAndRecParticleBase, ProtTomoBase):
    """Make pseudo-subtomograms"""

    _label = 'Make pseudo-subtomograms'
    _possibleOutputs = outputObjects

    def __int__(self, **kwargs):
        super().__int__(**kwargs)

    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        super()._defineParams(form)
        form.addSection(label='Reconstruct pseudo-Subtomograms')
        super()._defineCommonRecParams(form)
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
                      default=10,
                      help='It is the (full) opening angle of the cone to be suppressed, given in degrees. This angle '
                           'should  include both the uncertainty about the membrane orientation and its variation '
                           'across the region represented in the subtomogram.')
        form.addParam('applyOffsets', BooleanParam,
                      label='Apply offsets?',
                      default=False,
                      help='If set to Yes, rlnOrigin<X/Y/Z> translations are combined with rlnCoordinate<X/Y/Z> to '
                           'construct subtomos on their refined centers.')
        form.addParam('outputInFloat16', BooleanParam,
                      label='Write output in float16?',
                      default=True,
                      help='If set to Yes, this program will write output images in float16 MRC format. This will '
                           'save a factor of two in disk space compared to the default of writing in float32. Note '
                           'that RELION and CCPEM will read float16 images, but other programs may not (yet) do so.')

    # -------------------------- INSERT steps functions -----------------------
    def _insertAllSteps(self):
        self._insertFunctionStep(self.convertInputStep)
        self._insertFunctionStep(self.relionMakePseudoSubtomos)
        self._insertFunctionStep(self.createOutputStep)

    # -------------------------- STEPS functions ------------------------------
    def convertInputStep(self):
        writeSetOfPseudoSubtomograms(self.inReParticles.get(), self.getOutStarFile())

    def relionMakePseudoSubtomos(self):
        Plugin.runRelionTomo(self, getProgram('relion_tomo_subtomo', self.numberOfMpi.get()),
                             self._genMakePseudoSubtomoCmd(), numberOfMpi=self.numberOfMpi.get())

    def createOutputStep(self):
        inReParticles = self.inReParticles.get()
        # Pseudosubtomos
        psubtomoSet = createSetOfRelionPSubtomograms(self._getPath(),
                                                     self._getExtraPath(OPTIMISATION_SET_STAR),
                                                     template=PSUBTOMOS_SQLITE,
                                                     tsSamplingRate=inReParticles.getSamplingRate(),
                                                     relionBinning=self.binningFactor.get(),
                                                     boxSize=self.croppedBoxSize.get())
        # Fill the set with the generated particles
        readSetOfPseudoSubtomograms(psubtomoSet)

        self._defineOutputs(**{outputObjects.relionParticles.name: psubtomoSet})
        self._defineSourceRelation(inReParticles, psubtomoSet)

    # -------------------------- INFO functions -------------------------------

    # # --------------------------- UTILS functions -----------------------------
    def _genMakePseudoSubtomoCmd(self):
        cmd = self._genCommonCmd()
        cmd += '--o %s ' % self._getExtraPath()
        if self.applyConeWeight.get():
            cmd += '--cone_weight '
            cmd += '--cone_angle %.2f ' % self.coneAngle.get()
        if self.applyOffsets.get():
            cmd += '--apply_offsets '
        if self.outputInFloat16.get():
            cmd += '--float16 '
        return cmd
