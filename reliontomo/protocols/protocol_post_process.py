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
from pwem.objects import VolumeMask, FSC
from pyworkflow.object import String
from pyworkflow.protocol import PointerParam, BooleanParam, FloatParam, GE, LE, IntParam, FileParam
from pyworkflow.utils import makePath, Message
from reliontomo import Plugin
from reliontomo.constants import POST_PROCESS_MRC, POSTPROCESS_DIR, \
    POSTPROCESS_STAR_FSC_TABLE, \
    POSTPROCESS_STAR_FSC_COLUMNS, FSC_REF_STAR, POSTPROCESS_STAR_FIELD
from reliontomo.protocols.protocol_base_relion import ProtRelionTomoBase

NO_MTF_FILE = 0


class outputObjects(Enum):
    postProcessVolume = VolumeMask
    outputFSC = FSC


class ProtRelionPostProcess(ProtRelionTomoBase):
    """Sharpen a 3D reference map and estimate the gold-standard FSC curves for subtomogram averaging
    After performing a refinement, the map needs to be sharpened. Also, the gold-standard FSC curves
    inside the auto-refine procedures only use unmasked maps (unless you’ve used the option
    Use solvent-flattened FSCs). This means that the actual resolution is under-estimated during
    the actual refinement, because noise in the solvent region will lower the FSC curve. relion’s
    procedure for B-factor sharpening and calculating masked FSC curves [CMF+13] is called post-processing.
    First however, we’ll need to make a mask to define where the protein ends and the solvent region starts.
    """

    _label = 'Post-processing'
    _possibleOutputs = outputObjects

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        # super()._defineCommonInputParams(form)
        form.addSection(label=Message.LABEL_INPUT)
        form.addParam('inVolume', PointerParam,
                      pointerClass='Volume',
                      label='Volume to sharpen',
                      important=True,
                      help='It will provide the two unfiltered half-reconstructions that were output upon convergence '
                           'of a 3D auto-refine run.')
        form.addParam('solventMask', PointerParam,
                      pointerClass='VolumeMask',
                      important=True,
                      label='Solvent mask',
                      help='Provide a soft mask where the protein is white (1) and the solvent is black (0). '
                           'Often, the softer the mask the higher resolution estimates you will get. A soft '
                           'edge of 5-10 pixels is often a good edge width.')
        form.addParam('calPixSize', FloatParam,
                      default=-1,
                      label='Calibrated pixel size (Å/px)',
                      help='Provide the final, calibrated pixel size in Angstroms. This value may be different from '
                           'the pixel-size used thus far, e.g. when you have recalibrated the pixel size using the fit '
                           'to a PDB model. The X-axis of the output FSC plot will use this calibrated value.')
        form.addSection(label='Sharpen')
        form.addParam('estimateBFactor', BooleanParam,
                      label='Estimate B-factor automatically?',
                      default=True,
                      help='If set to Yes, then the program will use the automated procedure described by Rosenthal '
                           'and Henderson (2003, JMB) to estimate an overall B-factor for your map, and sharpen it '
                           'accordingly. Note that your map must extend well beyond the lowest resolution included '
                           'in the procedure below, which should not be set to resolutions much lower than 10 '
                           'Angstroms. ')
        form.addParam('lowestResBFit', FloatParam,
                      default=10,
                      validators=[GE(10), LE(15)],
                      condition='estimateBFactor',
                      label='Lowest resolution for auto-B factor',
                      help='This is the lowest frequency (in Angstroms) that will be included in the linear fit of '
                           'the Guinier plot as described in Rosenthal and Henderson (2003, JMB). Dont use values '
                           'much lower or higher than 10 Angstroms. If your map does not extend beyond 10 Angstroms, '
                           'then instead of the automated procedure use your own B-factor.')
        form.addParam('useOwnBFactor', BooleanParam,
                      default=False,
                      label='Use your own B-factor?',
                      help='Instead of using the automated B-factor estimation, provide your own value. Use negative '
                           'values for sharpening the map. This option is useful if your map does not extend beyond '
                           'the 10A needed for the automated procedure, or when the automated procedure does not give '
                           'a suitable value (e.g. in more disordered parts of the map).')
        form.addParam('userBFactor', IntParam,
                      default=-1000,
                      validators=[LE(0)],
                      label='User-provided B-factor',
                      condition='useOwnBFactor',
                      help='Use negative values for sharpening. Be careful: if you over-sharpen your map, you may '
                           'end up interpreting noise for signal!')
        form.addParam('skipFscWeight', BooleanParam,
                      default=False,
                      label='Skip FSC-weighting?',
                      help='If set to No, then the output map will be low-pass filtered according to the '
                           'mask-corrected, gold-standard FSC-curve. Sometimes, it is also useful to provide an '
                           'ad-hoc low-pass filter, as due to local resolution variations some parts of the map may '
                           'be better and other parts may be worse than the overall resolution as measured by the FSC. '
                           'In such cases, set this option to Yes and provide an ad-hoc filter as described below.')
        form.addParam('adHocLowPassFilter', IntParam,
                      default=5,
                      validators=[GE(1)],
                      label='Ad-hoc low-pass filter (Å)',
                      condition='skipFscWeight',
                      help='This option allows one to low-pass filter the map at a user-provided frequency (in '
                           'Angstroms). When using a resolution that is higher than the gold-standard FSC-reported '
                           'resolution, take care not to interpret noise in the map for signal.')
        form.addParam('mtf', FileParam,
                      label='MTF of the detector',
                      help='If you know the MTF of your detector, provide it here. Curves for some well-known detectors'
                           ' may be downloaded from the RELION Wiki \n'
                           '- [[https://www3.mrc-lmb.cam.ac.uk/relion/index.php/ \n'
                           'Also see there for the exact format of your detector.  If you do not know the MTF of your'
                           ' detector and do not want to measure it, then by leaving this entry empty, you include '
                           'the MTF of your detector in your overall estimated B-factor upon sharpening the map.'
                           'Although that is probably slightly less accurate, the overall quality of your map '
                           'will probably not suffer very much.')
        form.addParam('origDetectorPixSize', FloatParam,
                      default=1,
                      validators=[GE(0.1), LE(2)],
                      condition='mtf',
                      label='Original detector pixel size ((Å)/pix)',
                      help='This is the original pixel size (in Angstroms) in the raw (non-super-resolution!) '
                           'micrographs.')
        self._defineExtraParams(form)

    # -------------------------- INSERT steps functions -----------------------
    def _insertAllSteps(self):
        makePath(self._getExtraPath(POSTPROCESS_DIR))
        self._insertFunctionStep(self.relionPostProcessStep)
        self._insertFunctionStep(self.createOutputStep)

    def relionPostProcessStep(self):
        Plugin.runRelionTomo(self, 'relion_postprocess', self.genPostProcessCmd())

    def createOutputStep(self):
        inVolume = self.inVolume
        fn = self._getExtraPath(FSC_REF_STAR)
        postProccesMrc = self._genPostProcessOutputMrcFile(POST_PROCESS_MRC)
        # Extend the sharpened volume with an attribute containing the postprocess.star file
        setattr(postProccesMrc, POSTPROCESS_STAR_FIELD, String(fn))

        # Output FSC
        setOfFSC = self.genFSCs(fn, POSTPROCESS_STAR_FSC_TABLE,
                                POSTPROCESS_STAR_FSC_COLUMNS)

        self._defineOutputs(**{outputObjects.postProcessVolume.name: postProccesMrc,
                               outputObjects.outputFSC.name: setOfFSC})
        self._defineSourceRelation(inVolume, postProccesMrc)

    # -------------------------- INFO functions -------------------------------
    def _validate(self):
        errorMsg = []
        if not self.inVolume.get().getHalfMaps():
            errorMsg.append('The introduced volume needs to get the corresponding half maps.')
        return errorMsg

    # --------------------------- UTILS functions -----------------------------
    def genPostProcessCmd(self):
        calPixSize = self.calPixSize.get() if self.calPixSize.get() > 0 else self.inVolume.get().getSamplingRate()
        half1, half2 = self.inVolume.get().getHalfMaps().split(',')
        cmd = ''
        cmd += '--i %s ' % half1
        cmd += '--i2 %s ' % half2
        cmd += '--o %s ' % self._getExtraPath(POSTPROCESS_DIR, POSTPROCESS_DIR.lower())
        cmd += '--mask %s ' % self.solventMask.get().getFileName()
        cmd += '--angpix %.2f ' % calPixSize
        # Sharpening
        if self.mtf.get():
            cmd += '--mtf %s ' % self.mtf.get()
            cmd += '--mtf_angpix %.2f' % self.origDetectorPixSize.get()
        if self.estimateBFactor.get():
            cmd += '--auto_bfac --autob_lowres %.2f ' % self.lowestResBFit.get()
        if self.useOwnBFactor.get():
            cmd += '--adhoc_bfac %.2f ' % self.userBFactor.get()
        # Filtering
        if self.skipFscWeight.get():
            cmd += '--skip_fsc_weighting --low_pass %i ' % self.adHocLowPassFilter.get()
        # Extra params
        cmd += self._genExtraParamsCmd()

        return cmd
