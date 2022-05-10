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

from pwem.objects import VolumeMask
from pwem.protocols import EMProtocol
from pyworkflow import BETA
from pyworkflow.protocol import PointerParam, BooleanParam, FloatParam, GE, LE, IntParam, FileParam
from pyworkflow.utils import Message, makePath
from reliontomo import Plugin
from reliontomo.constants import POST_PROCESS_MRC, POST_PROCESS_MASKED_MRC, POSTPROCESS_DIR
from reliontomo.objects import relionTomoMetadata
from reliontomo.utils import genRelionParticles

NO_MTF_FILE = 0


class outputObjects(Enum):
    relionParticles = relionTomoMetadata
    postProcessVolume = VolumeMask


class ProtRelionPostProcess(EMProtocol):
    """Sharpen a 3D reference map and estimate the gold-standard FSC curves for subtomogram averaging"""

    _label = 'Sharpen a 3D reference maps'
    _devStatus = BETA
    _possibleOutputs = outputObjects

    def __init__(self, **kargs):
        super().__init__(**kargs)

    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        form.addSection(label=Message.LABEL_INPUT)
        form.addParam('inOptSet', PointerParam,
                      pointerClass='relionTomoMetadata',
                      label='Input Relion Tomo Metadata',
                      important=True)
        form.addParam('inVolume', PointerParam,
                      pointerClass='AverageSubTomogram',
                      label='Use halves from this refined volume',
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
                      default=3,
                      validators=[GE(0.3), LE(5)],
                      label='Calibrated pixel size (Å/pix)',
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
                      label='Lowest resolution for auto-B fit',
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
                      validators=[GE(-2000), LE(0)],
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
                      validators=[GE(1), LE(40)],
                      label='Ad-hoc low-pass filter (Å)',
                      condition='skipFscWeight',
                      help='This option allows one to low-pass filter the map at a user-provided frequency (in '
                           'Angstroms). When using a resolution that is higher than the gold-standard FSC-reported '
                           'resolution, take care not to interpret noise in the map for signal.')
        form.addParam('mtf', FileParam,
                      label='MTF of the detector',
                      help='User-provided STAR-file with the MTF-curve '
                           'of the detector. Use the wizard to load one '
                           'of the predefined ones provided at:\n'
                           '- [[https://www3.mrc-lmb.cam.ac.uk/relion/index.php/'
                           'FAQs#Where_can_I_find_MTF_curves_for_typical_detectors.3F]'
                           '[Relion\'s Wiki FAQs]]\n'
                           ' - [[https://www.gatan.com/techniques/cryo-em#MTF][Gatan\'s website]]\n\n'
                           'Relion param: *--mtf*')
        form.addParam('origDetectorPixSize', FloatParam,
                      default=1,
                      validators=[GE(0.1), LE(2)],
                      condition='mtf',
                      label='Original detector pixel size ((Å)/pix)',
                      help='This is the original pixel size (in Angstroms) in the raw (non-super-resolution!) '
                           'micrographs.')

    # -------------------------- INSERT steps functions -----------------------
    def _insertAllSteps(self):
        makePath(self._getExtraPath(POSTPROCESS_DIR))
        self._insertFunctionStep(self.relionPostProcessStep)
        self._insertFunctionStep(self.createOutputStep)

    def relionPostProcessStep(self):
        Plugin.runRelionTomo(self, 'relion_postprocess', self.genPostProcessCmd())

    def createOutputStep(self):
        inOptSet = self.inOptSet.get()
        # Output RelionParticles
        relionParticles = genRelionParticles(self._getExtraPath(), inOptSet)
        # Output FSC masks
        postProccesMrc = self._genPostProcessOutputMrcFile(POST_PROCESS_MRC)
        # postProcessMaskedMrc = self._genPostProcessOutputMrcFile(POST_PROCESS_MASKED_MRC)

        outputDict = {outputObjects.relionParticles.name: relionParticles,
                      outputObjects.postProcessVolume.name: postProccesMrc}
        self._defineOutputs(**outputDict)
        self._defineSourceRelation(inOptSet, relionParticles)
        self._defineSourceRelation(inOptSet, postProccesMrc)

    # -------------------------- INFO functions -------------------------------
    def _validate(self):
        pass

    # --------------------------- UTILS functions -----------------------------
    def genPostProcessCmd(self):
        half1, half2 = self.inVolume.get().getHalfMaps().split(',')
        cmd = ''
        cmd += '--i %s ' % half1
        cmd += '--i2 %s ' % half2
        cmd += '--o %s ' % self._getExtraPath(POSTPROCESS_DIR, POSTPROCESS_DIR.lower())
        cmd += '--mask %s ' % self.solventMask.get().getFileName()
        cmd += '--angpix %.2f ' % self.calPixSize.get()
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

        return cmd

    def _genPostProcessOutputMrcFile(self, fileName):
        postProccesMrc = VolumeMask()
        postProccesMrc.setFileName(self._getExtraPath(POSTPROCESS_DIR, fileName))
        postProccesMrc.setSamplingRate(self.inOptSet.get().getCurrentSamplingRate())

        return postProccesMrc
