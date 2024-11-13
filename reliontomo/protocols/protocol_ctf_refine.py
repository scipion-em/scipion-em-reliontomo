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
from pyworkflow.protocol import IntParam, BooleanParam, GE, LE, FloatParam, EnumParam
from reliontomo import Plugin
from reliontomo.protocols.protocol_base_per_part_per_tilt import ProtRelionPerParticlePerTiltBase
from reliontomo.protocols.protocol_base_relion import IS_RELION_50
from reliontomo.utils import getProgram

oddAberrationOrders = [3, 5, 7]
evenAberrationOrders = [4, 6, 8]


class ProtRelionCtfRefine(ProtRelionPerParticlePerTiltBase):
    """Tomo CTF refine:\n

    This program estimates the astigmatic defoci of the individual tilt images, the ice
    thickness (which determines the overall signal intensity) and higher-order optical aberrations.\n

    _Defocus_: In tomography, the relative depth distances between particles are known from the 3D
    positions of the particles. Therefore, only one defocus value is estimated for all the particles
    in each tilt image. Because of the often large number of particles in each tomogram, this value
    can typically be estimated to greater precision than in single-particle analysis, where the defocus
    of each particle has to be estimated independently.\n
    _Ice-thickness_: A thicker sample permits fewer electrons to pass through, which reduces the scale
    of the signal. In addition to the actual variation in sample thickness, the effective thickness
    of the ice also increases as the sample is tilted. This program allows the user to estimate the
    signal intensity either independently for each tilt image, or by fitting the base thickness,
    the initial beam luminance and the surface normal of the sample assuming Beer-Lambert’s Law.
    In the latter case, only one surface normal and base thickness are estimated for an entire
    tilt series, which allows for a more stable fit from tomograms with fewer particles.\n
    _Higher-order optical aberrations_: This algorithm works analogously to relion_ctf_refine for
    single-particle analysis. As in single-particle analysis, the aberrations are estimated
    per optics group. This allows the user to group particles that are expected to share the
    same aberrations, either by image region or by subset of tilt series, or both. Both
    symmetrical (even) and antisymmetrical (odd) aberrations are supported. A detailed description
     of the aberrations estimation algorithm can be found in the relion aberrations paper.
    """

    _label = 'CTF refinement'

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        super()._defineParams(form)
        form.addSection(label='Defocus')
        super()._insertBoxSizeForEstimationParam(form)
        form.addParam('refineDefocus', BooleanParam,
                      label='Refine defocus?',
                      default=True,
                      help='If set to Yes, then estimate the defoci of the individual tilt images.')
        form.addParam('defocusRange', IntParam,
                      label="Defocus search range (Å)",
                      condition='refineDefocus',
                      default=3000,
                      validators=[GE(0)],
                      help='Defocus search range (in A). This search range will be, by default, '
                           'sampled in 100 steps. Use the additional argument --ds to change the '
                           'number of sampling points.')
        form.addParam('doDefocusReg', BooleanParam,
                      label="Do defocus regularisation?",
                      condition='refineDefocus',
                      default=False,
                      help="Apply defocus regularisation.\n\nHigh-tilt images do not offer enough signal to recover "
                           "the defocus value precisely. The regularisation forces the estimated defoci to assume "
                           "similar values within a given tilt series, which prevents those high-tilt images from "
                           "overfitting.")
        form.addParam('regParam', FloatParam,
                      label="Defocus regularisation scale",
                      condition='doDefocusReg',
                      default=0.1,
                      validators=[GE(0), LE(1)],
                      help='This is the strength of the defocus regularizer')
        form.addParam('refineContrast', BooleanParam,
                      label='Refine contrast scale?',
                      default=True,
                      help='If set to Yes, then estimate the signal scale or ice thickness.')
        form.addParam('refineScalePerFrame', BooleanParam,
                      label='Refine scale per frame?',
                      default=True,
                      help='If set to Yes, then estimate the signal-scale parameter independently for each tilt. If '
                           'not specified, the ice thickness, beam luminance and surface normal are estimated instead. '
                           'Those three parameters then imply the signal intensity for each frame. Due to the smaller '
                           'number of parameters, the ice thickness model is more robust to noise. By default, the ice '
                           'thickness and surface normal will be estimated per tilt-series, and the beam luminance '
                           'globally.')
        form.addParam('refineScalePerTomo', BooleanParam,
                      label="Refine scale per tomogram?",
                      default=False,
                      help="If set to Yes, then the beam luminance will be estimated separately for each tilt series. "
                           "This is not recommended.")
        if not IS_RELION_50:
            form.addSection(label='Aberrations')
            form.addParam('refineOddAbe', BooleanParam,
                          label="Refine odd aberrations?",
                          default=True,
                          help="If set to Yes, then odd higher-order aberrations will be estimated. These are"
                               "the asymmetrical aberrations")
            form.addParam('maxAbeOddOrder', EnumParam,
                          display=EnumParam.DISPLAY_HLIST,
                          label='Max order of odd aberrations',
                          condition='refineOddAbe',
                          choices=oddAberrationOrders,
                          default=0,
                          help='The third order aberration will be comma and trefoil. Higer aberrations as pentafoil '
                               'are barely considered')
            form.addParam('refineEvenAbe', BooleanParam,
                          label="Refine even aberrations?",
                          default=True,
                          help="If set to Yes, then even higher-order aberrations will be estimated. These are"
                               "the symmetrical aberrations")
            form.addParam('maxAbeEvenOrder', EnumParam,
                          display=EnumParam.DISPLAY_HLIST,
                          label='Max order of even aberrations',
                          condition='refineEvenAbe',
                          choices=evenAberrationOrders,
                          default=0,
                          help='The forth order aberrations are spherical aberration, quadrafoil and secondary '
                               'astigmatism, higher aberrations are barely considered.')

        self._defineExtraParams(form)
        form.addParallelSection(threads=1, mpi=1)

    # -------------------------- INSERT steps functions -----------------------
    def _insertAllSteps(self):
        self._insertFunctionStep(self.convertInputStep, needsGPU=False)
        self._insertFunctionStep(self._relionCTFRefine, needsGPU=False)
        self._insertFunctionStep(self.createOutputStep, needsGPU=False)

    # -------------------------- STEPS functions ------------------------------
    def _relionCTFRefine(self):
        Plugin.runRelionTomo(self, getProgram('relion_tomo_refine_ctf', self.numberOfMpi.get()),
                             self._genTomoRefineCtfCmd(), numberOfMpi=self.numberOfMpi.get())

    # -------------------------- INFO functions -------------------------------
    def _validate(self):
        errorMsgs = []
        if self.refineContrast.get():
            if self.refineScalePerFrame.get() and self.refineScalePerTomo.get():
                errorMsgs.append('Per-tomogram scale estimation and per-frame scale estimation are mutually exclusive')
        return errorMsgs

    # --------------------------- UTILS functions -----------------------------
    def _genTomoRefineCtfCmd(self):
        cmd = self._genIOCommand()
        if self.refineDefocus.get():
            cmd += '--do_defocus '
            cmd += '--d0 %i ' % self.defocusRange.get()
            cmd += '--d1 %i ' % self.defocusRange.get()
            if self.doDefocusReg.get():
                cmd += '--do_reg_defocus --lambda %.2f ' % self.regParam.get()
        if self.refineContrast.get():
            cmd += '--do_scale '
            if self.refineScalePerFrame.get():
                cmd += '--per_frame_scale '
            if self.refineScalePerTomo.get():
                cmd += '--per_tomo_scale '
        if not IS_RELION_50:
            if self.refineEvenAbe.get():
                cmd += '--do_even_aberrations --ne %i ' % evenAberrationOrders[self.maxAbeEvenOrder.get()]
            if self.refineOddAbe.get():
                cmd += '--do_odd_aberrations --no %i ' % oddAberrationOrders[self.maxAbeOddOrder.get()]

        return cmd
