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
    """
    Refines tomographic Contrast Transfer Function parameters from tilt-series data in order to improve
    the accuracy and interpretability of subtomogram reconstructions. The protocol estimates defocus,
    signal attenuation effects related to ice thickness, and higher-order optical aberrations across
    tilt images, allowing the resulting tomographic data to achieve improved contrast consistency and
    optical correction.

    AI Generated:

    Tomo CTF Refine (ProtRelionCtfRefine) — User Manual
        Overview

        The Tomo CTF Refine protocol performs advanced refinement of optical parameters in cryo-electron
        tomography datasets. Its purpose is to improve the physical description of image formation within
        a tilt series by refining parameters that directly influence contrast transfer and signal quality.
        Accurate correction of these effects is essential for achieving higher-resolution subtomogram
        averages and for improving the reliability of downstream structural interpretation.

        In tomography workflows, particles occupy different depths inside the specimen and are observed
        across many tilted projections. Because of this geometry, optical distortions and signal
        attenuation vary systematically throughout the dataset. This protocol takes advantage of the
        known three-dimensional particle positions to estimate more stable and physically meaningful
        corrections than those typically achievable in single-particle analysis.

        Defocus Refinement

        One of the primary functions of the protocol is the refinement of defocus values for tilt
        images. In cryo-electron tomography, all particles from the same projection share a common
        imaging geometry, which allows defocus estimation to be performed collectively for each tilt.
        Since many particles contribute simultaneously, the resulting estimates are often more precise
        and robust than independent particle-by-particle approaches.

        Defocus refinement is especially important for high-resolution subtomogram averaging, where
        small inaccuracies in phase correction can significantly reduce map quality. The protocol allows
        users to define a search range for the refinement, enabling adaptation to datasets with varying
        imaging conditions or uncertainty in the initial microscope parameters.

        For challenging datasets acquired at high tilt angles, the protocol also supports regularization
        of defocus estimates. This stabilizes the refinement by constraining neighboring tilt images to
        maintain physically consistent values. Such regularization is particularly useful when the
        signal-to-noise ratio decreases strongly at high tilts, where overfitting can otherwise occur.

        Ice Thickness and Signal Scale Estimation

        The protocol also estimates signal attenuation effects caused by specimen thickness and electron
        absorption. In cryo-electron tomography, thicker ice regions reduce transmitted electron signal,
        leading to variations in image intensity across the tilt series. These effects become more
        pronounced at high tilt angles because electrons traverse a longer effective path through the
        sample.

        Users may choose between estimating signal intensity independently for each tilt image or using
        a physically constrained ice-thickness model. Independent estimation provides flexibility and
        can adapt to complex datasets with irregular signal variations. In contrast, the ice-thickness
        model estimates global physical parameters such as specimen thickness, beam luminance, and
        surface orientation, resulting in a more stable and noise-resistant refinement.

        The physically constrained model is often preferable for datasets with limited particle numbers
        or noisy projections because it reduces the number of free parameters. However, datasets with
        strong local variations in illumination or sample thickness may benefit from the more flexible
        per-frame refinement approach.

        Higher-Order Optical Aberration Refinement

        The protocol optionally refines higher-order optical aberrations that arise from imperfections
        in the microscope optics. These aberrations include both symmetrical and asymmetrical distortions
        that affect image quality and limit achievable resolution if left uncorrected.

        Symmetrical aberrations include effects such as spherical aberration and related higher-order
        distortions, while asymmetrical aberrations include coma-like or trefoil-like components. The
        protocol estimates these parameters per optics group, allowing users to separate datasets
        acquired under different optical conditions or from different microscope regions.

        In practice, higher-order aberration refinement is most beneficial for high-resolution datasets
        where subtle optical imperfections become resolution limiting. For lower-resolution exploratory
        projects, these corrections may provide smaller improvements and can often be omitted during
        initial processing stages.

        Inputs and General Workflow

        The protocol requires aligned tomographic particle data together with previously estimated CTF
        information. Reliable particle alignment and accurate tilt-series geometry are essential because
        the refinement depends directly on the spatial consistency of the dataset.

        During execution, the protocol analyzes the tilt-series projections collectively and estimates
        the selected optical parameters according to the chosen refinement strategy. Depending on the
        dataset size and selected options, the refinement may range from relatively fast defocus updates
        to more computationally demanding aberration estimation procedures.

        The resulting refined optical parameters can then be propagated into subsequent subtomogram
        averaging or refinement workflows, leading to improved contrast correction and potentially
        higher-resolution reconstructions.

        Practical Recommendations

        For most biological projects, beginning with defocus refinement and signal-scale estimation is
        a sensible strategy. These corrections often provide the largest immediate improvements in map
        quality and alignment stability. Defocus regularization is strongly recommended for datasets
        containing many high-tilt images or low signal levels.

        Ice-thickness modeling is particularly useful for datasets with limited particle counts because
        it imposes physically meaningful constraints that improve robustness. Per-frame scale refinement
        may be preferable when illumination conditions vary substantially across the acquisition.

        Higher-order aberration refinement is generally most valuable during late-stage refinement of
        high-resolution datasets. Users should ensure that optics groups are defined appropriately so
        that particles sharing similar optical conditions are refined together.

        Final Perspective

        Accurate modeling of optical effects is one of the key factors determining the quality of
        cryo-electron tomography reconstructions. By refining defocus, correcting signal attenuation,
        and compensating for higher-order aberrations, this protocol improves the physical consistency
        of tomographic datasets and enhances the reliability of downstream structural interpretation.
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
        form.addParallelSection(threads=0, mpi=1)

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
