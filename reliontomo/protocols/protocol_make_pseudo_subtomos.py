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
from reliontomo import Plugin
from reliontomo.constants import OUT_SUBTOMOS_STAR
from reliontomo.convert import readSetOfPseudoSubtomograms
from reliontomo.objects import SetOfPseudoSubtomograms
from reliontomo.protocols.protocol_base_make_pseusosubtomos_and_rec_particle import \
    ProtRelionMakePseudoSubtomoAndRecParticleBase
from tomo.protocols import ProtTomoBase


class ProtRelionMakePseudoSubtomograms(ProtRelionMakePseudoSubtomoAndRecParticleBase, ProtTomoBase):
    """Make pseudo-subtomograms"""

    _label = 'Make pseudo-subtomograms'

    def __init__(self, **args):
        ProtRelionMakePseudoSubtomoAndRecParticleBase.__init__(self, **args)

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
        form.addParam('coneWidthZ0', IntParam,
                      label='Cone width at Z = 0 (pix.)',
                      condition='applyConeWeight',
                      default=2,
                      expertLevel=LEVEL_ADVANCED)
        form.addParam('keepParticlesDark', BooleanParam,
                      label='Do not invert contrast (keep particles dark)',
                      default=False,
                      expertLevel=LEVEL_ADVANCED)
        form.addParam('write3DCtfs', BooleanParam,
                      label='Write 3D CTFs?',
                      default=False,
                      expertLevel=LEVEL_ADVANCED)
        form.addParam('writeCtfCorrectedSubtomos', BooleanParam,
                      label='Write CTF-corrected subtomograms?',
                      default=False,
                      expertLevel=LEVEL_ADVANCED)

    # -------------------------- INSERT steps functions -----------------------
    def _insertAllSteps(self):
        self._insertFunctionStep(self._relionMakePseudoSubtomos)
        self._insertFunctionStep(self.createOutputStep)

    # -------------------------- STEPS functions ------------------------------
    def _relionMakePseudoSubtomos(self):
        Plugin.runRelionTomo(self, 'relion_tomo_subtomo', self._genMakePseudoSubtomoCmd(),
                             numberOfMpi=self.numberOfMpi.get())

    def createOutputStep(self):
        starFile = self._getExtraPath(OUT_SUBTOMOS_STAR)
        protPrepDataInputSubtomos = self.inputPrepareDataProt.get().inputSubtomos.get()

        outputSet = self._createSet(SetOfPseudoSubtomograms, 'pseudosubtomograms%s.sqlite', '')
        outputSet.setOpticsGroupObjFromStar(starFile)
        outputSet.setSamplingRate(protPrepDataInputSubtomos.getSamplingRate())
        precedents = protPrepDataInputSubtomos.getCoordinates3D().get().getPrecedents()
        readSetOfPseudoSubtomograms(starFile, precedents, outputSet, invert=True)
        self._defineOutputs(outputSetOfPseudoSubtomogram=outputSet)

    # # -------------------------- INFO functions -------------------------------
    def _validate(self):
        pass

    # # --------------------------- UTILS functions -----------------------------
    def _genMakePseudoSubtomoCmd(self):
        cmd = self._genCommonCmd()
        cmd += '--o %s ' % self._getExtraPath()
        if self.applyConeWeight.get():
            cmd += '--cone_weight '
            cmd += '--cone_angle %.2f ' % self.coneAngle.get()
            cmd += '--cone_sig0 %i ' % self.coneWidthZ0.get()
        if self.keepParticlesDark.get():
            cmd += ' --no_ic '
        if self.write3DCtfs.get():
            cmd += '--ctf '
        if self.writeCtfCorrectedSubtomos.get():
            cmd += '--div '
        return cmd
