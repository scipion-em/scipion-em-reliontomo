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
from pyworkflow.protocol import LEVEL_ADVANCED, IntParam, FloatParam, BooleanParam
from relion.convert import OpticsGroups
from reliontomo import Plugin
from reliontomo.constants import OUT_SUBTOMOS_STAR
from reliontomo.convert import readSetOfPseudoSubtomograms
from reliontomo.protocols.protocol_base_make_pseusosubtomos_and_rec_particle import \
    ProtRelionMakePseudoSubtomoAndRecParticleBase
from reliontomo.utils import getProgram
from tomo.objects import SetOfSubTomograms
from tomo.protocols import ProtTomoBase


class ProtRelionMakePseudoSubtomograms(ProtRelionMakePseudoSubtomoAndRecParticleBase, ProtTomoBase):
    """Make pseudo-subtomograms"""

    _label = 'Make pseudo-subtomograms'

    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        ProtRelionMakePseudoSubtomoAndRecParticleBase._defineParams(self, form)
        form.addSection(label='Reconstruct pseudo-Subtomograms')
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
        form.addParam('rotateVolumes', BooleanParam,
                      label='Rotate volumes?',
                      default=False,
                      help='If set to Yes, rlnAngle<Rot/Tilt/Psi> orientations are added to '
                           'rlnTomoSubtomogram<Rot/Tilt/Psi> to construct rotated volumes.')
        form.addParam('restoreOrientations', BooleanParam,
                      label='Restore orientations?',
                      default=False,
                      help='If set to Yes, rlnTomoSubtomogram<Rot/Tilt/Psi> orientations are added to '
                           'rlnAngle<Rot/Tilt/Psi>. Particles are not constructed.')
        form.addParam('outputInFloat16', BooleanParam,
                      label='Write output in float16?',
                      default=True,
                      help='If set to Yes, this program will write output images in float16 MRC format. This will '
                           'save a factor of two in disk space compared to the default of writing in float32. Note '
                           'that RELION and CCPEM will read float16 images, but other programs may not (yet) do so.')

    # -------------------------- INSERT steps functions -----------------------
    def _insertAllSteps(self):
        self._initialize()
        self._insertFunctionStep(self.convertInputStep)
        self._insertFunctionStep(self._relionMakePseudoSubtomos)
        self._insertFunctionStep(self.createOutputStep)

    # -------------------------- STEPS functions ------------------------------

    def _relionMakePseudoSubtomos(self):
        Plugin.runRelionTomo(self, getProgram('relion_tomo_subtomo', self.numberOfMpi.get()),
                             self._genMakePseudoSubtomoCmd(), numberOfMpi=self.numberOfMpi.get())

    def createOutputStep(self):
        starFile = self._getExtraPath(OUT_SUBTOMOS_STAR)
        protPrepData = self.inputPrepareDataProt.get()
        tsSamplingRate = protPrepData.inputCtfTs.get().getSetOfTiltSeries().getSamplingRate()  # bin1
        outputSet = self._createSet(SetOfSubTomograms, 'subtomograms%s.sqlite', '')
        outputSet.getAcquisition().opticsGroupInfo.set(OpticsGroups.fromStar(starFile).toString())
        outputSet.setSamplingRate(tsSamplingRate * self.binningFactor.get())
        readSetOfPseudoSubtomograms(starFile, outputSet)
        self._defineOutputs(outputPseudoSubtomograms=outputSet)

    # # -------------------------- INFO functions -------------------------------
    def _validate(self):
        validationMsg = []
        if self.rotateVolumes.get() and self.restoreOrientations.get():
            validationMsg.append('Restore angles and rotate p-subtomos cannot be applied simultaneously.')

        return validationMsg

    # # --------------------------- UTILS functions -----------------------------
    def _genMakePseudoSubtomoCmd(self):
        cmd = self._genCommonCmd()
        cmd += '--o %s ' % self._getExtraPath()
        if self.applyConeWeight.get():
            cmd += '--cone_weight '
            cmd += '--cone_angle %.2f ' % self.coneAngle.get()
        if self.outputInFloat16.get():
            cmd += '--float16 '
        if self.rotateVolumes.get():
            cmd += '--apply_angles '
        if self.restoreOrientations.get():
            cmd += '--restore '
        return cmd
