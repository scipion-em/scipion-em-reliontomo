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
    """Make pseudo-subtomograms\n

    Pseudo-subtomograms do not aim to accurately represent the scattering potential of the underlying particles.
    Instead, they serve as a practical means to implement an approximation to the 2D approach within
    the existing RELION framework. In the original RELION4 article an accurate defition is given, see: \n

     https://doi.org/10.7554/eLife.83724
    \n
    Before processing with any of the regular relion programs (those not specifically intended for
    tomography, e.g. relion_refine), we first need to construct the individual pseudo-subtomogram
    particles, equivalent to the particle extraction process in the SPA workflow.\

    A more technical explanation, pseudo-subtomograms are 3D-Arrays (volumes) that constructed from
    the sums of 2D tilt-series images pre-multiplied by contrast transfer functions (CTFs),
    along with auxiliary arrays that store the corresponding sum of squared CTFs and the frequency
    of observation for each 3D voxel.
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
