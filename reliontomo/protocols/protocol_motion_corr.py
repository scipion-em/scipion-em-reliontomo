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
from os.path import join

from pyworkflow.protocol import PointerParam, IntParam, GE, BooleanParam, LEVEL_ADVANCED, FloatParam, EnumParam, \
    FileParam
from pyworkflow.utils import Message, makePath
from reliontomo import Plugin
from reliontomo.constants import IN_TS_STAR
from reliontomo.convert import convert50_tomo
from reliontomo.protocols.protocol_base_relion import ProtRelionTomoBase
from reliontomo.utils import getProgram
from tomo.objects import SetOfTiltSeries

# Gain rotation values
NO_ROTATION = 0
DEG_90 = 1
DEG_180 = 2
DEF_270 = 3

# Gain flip values
NO_FLIP = 0
FLIP_UPSIDE_DOWN = 1
FLIP_LEFT_RIGHT = 2


class outputObjects(Enum):
    tiltSeries = SetOfTiltSeries


class ProtRelionTomoMotionCorr(ProtRelionTomoBase):
    """Motion correction"""
    _label = 'Motion correction'
    _possibleOutputs = outputObjects

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.tsDict = None

    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        form.addSection(label=Message.LABEL_INPUT)
        form.addParam('inputTiltSeriesM', PointerParam,
                      pointerClass='SetOfTiltSeriesM',
                      important=True,
                      label='Tilt-Series movies')
        form.addParam('outputInFloat16', BooleanParam,
                      label='Write output in float16?',
                      expertLevel=LEVEL_ADVANCED,
                      default=False,
                      help="If set to Yes, RelionCor2 will write output images in float16 MRC format. This will save a "
                           "factor of two in disk space compared to the default of writing in float32. Note that "
                           "RELION and CCPEM will read float16 images, but other programs may not (yet) do so. "
                           "Also note that this option does not work with UCSF MotionCor2.")
        form.addParam('saveEvenOdd', BooleanParam,
                      label='Save even/odd images?',
                      default=False,
                      help="It will write output images summed from both the even frames of the input movie and the "
                           "odd frames of the input movie. This generates two versions of the same movie which are "
                           "essential if you wish to carry out denoising later with cryoCARE. If you are unsure "
                           "whether you will need denoising later, it is best to select Yes, but be aware this option "
                           "increases the processing time.")
        form.addParam('saveSumPowerSpectra', BooleanParam,
                      label='Save sum of power spectra?',
                      default=True,
                      help="Sum of non-dose weighted power spectra provides better signal for CTF estimation. You "
                           "must use this option when writing in float16.")
        form.addParam('sumPowerSpectraNFrames', IntParam,
                      label='Sum power spectra every n frames',
                      default=4,
                      condition='saveSumPowerSpectra',
                      validators=[GE(0)],
                      help="McMullan et al (Ultramicroscopy, 2015) suggest summing power spectra every 4.0 e/A2 gives "
                           "optimal Thon rings.")

        form.addSection(label='Motion')
        form.addParam('bFactor', IntParam,
                      label='Bfactor',
                      default=150,
                      validators=[GE(0)],
                      help="The B-factor that will be applied to the micrographs.")
        line = form.addLine('Number of patches',
                            help='Number of patches to be used for patch based alignment.')
        line.addParam('patchX', IntParam, default=1, label='X')
        line.addParam('patchY', IntParam, default=1, label='Y')
        form.addParam('groupFrames', IntParam,
                      label='Group frames',
                      default=1,
                      validators=[GE(1)],
                      help="Average together this many frames before calculating the beam-induced shifts.")
        form.addParam('binningFactor', FloatParam,
                      label='Binning factor',
                      default=1,
                      validators=[GE(1)],
                      help="Bin the micrographs this much by a windowing operation in the Fourier Transform. Binning "
                           "at this level is hard to un-do later on, but may be useful to down-scale super-resolution "
                           "images. Float-values may be used. Do make sure though that the resulting micrograph size "
                           "is even.")
        form.addParam('gainRot', EnumParam,
                      default=NO_ROTATION,
                      choices=['No rotation (0)',
                               ' 90 degrees (1)',
                               '180 degrees (2)',
                               '270 degrees (3)'],
                      label='Gain rotation',
                      help="Rotate the gain reference by this number times 90 "
                           "degrees clockwise in relion_display. This is the "
                           "same as -RotGain in MotionCor2. \n"
                           "Note that MotionCor2 uses a different convention "
                           "for rotation so it says 'counter-clockwise'.")

        form.addParam('gainFlip', EnumParam,
                      default=NO_FLIP,
                      choices=['No flipping        (0)',
                               'Flip upside down   (1)',
                               'Flip left to right (2)'],
                      label='Gain flip',
                      help="Flip the gain reference after rotation. "
                           "This is the same as -FlipGain in MotionCor2. "
                           "0 means do nothing, 1 means flip Y (upside down) "
                           "and 2 means flip X (left to right).")

        form.addParam('defectFile', FileParam,
                      allowsNull=True,
                      label='Defects file',
                      help='Location of a UCSF MotionCor2-style '
                           'defect text file or a defect map that '
                           'describe the defect pixels on the detector. '
                           'Each line of a defect text file should contain '
                           'four numbers specifying x, y, width and height '
                           'of a defect region. A defect map is an image '
                           '(MRC or TIFF), where 0 means good and 1 means '
                           'bad pixels. The coordinate system is the same '
                           'as the input movie before application of '
                           'binning, rotation and/or flipping.\n\n'
                           '_Note that the format of the defect text is '
                           'DIFFERENT from the defect text produced '
                           'by SerialEM!_\n One can convert a SerialEM-style '
                           'defect file into a defect map using IMOD '
                           'utilities e.g.:\n'
                           '*clip defect -D defect.txt -f tif movie.tif defect_map.tif*\n'
                           'See explanations in the SerialEM manual.\n'
                           'Leave empty if you do not have any defects, '
                           'or do not want to correct for defects on your detector.')

        form.addSection(label='EER')
        form.addParam('eerFractionation', IntParam,
                      default=32,
                      validators=[GE(1)],
                      label='EER fractionation',
                      help="The number of hardware frames to group into one fraction. This option is relevant only "
                           "for Falcon4 movies in the EER format. Note that all 'frames' in the GUI (e.g. first and "
                           "last frame for corrected sum, dose per frame) refer to fractions, not raw detector frames. "
                           "See https://www3.mrc-lmb.cam.ac.uk/relion/index.php/Image_compression#Falcon4_EER for "
                           "detailed guidance on EER processing.")

        self._defineExtraParams(form)
        form.addParallelSection(threads=4, mpi=1)

    # -------------------------- INSERT steps functions -----------------------
    def _insertAllSteps(self):
        self._initialize()
        self._insertFunctionStep(self.convertInputStep)
        self._insertFunctionStep(self.correctMotionStep)
        # self._insertFunctionStep(self.createOutputStep)

    # -------------------------- STEPS functions ------------------------------
    def _initialize(self):
        inTsMSet = self.inputTiltSeriesM.get()
        presentTsIds = inTsMSet.getTSIds()
        self.tsDict = {tsM.getTsId(): tsM.clone(ignoreAttrs=[]) for tsM in inTsMSet if tsM.getTsId() in presentTsIds}

    def convertInputStep(self):
        inFilesPath = self.getInFilesPath()
        makePath(inFilesPath)
        writer = convert50_tomo.Writer()
        writer.tsMSet2Star(self.inputTiltSeriesM.get(), self.tsDict, inFilesPath)

    def correctMotionStep(self):
        nMpi = self.numberOfMpi.get()
        Plugin.runRelionTomo(self,
                             getProgram('relion_run_motioncorr', nMpi=nMpi),
                             self.getMotionCorrSubtomosCmd(),
                             numberOfMpi=nMpi)

    def createOutputStep(self):
        pass

    # -------------------------- INFO functions -------------------------------
    def _validate(self):
        errorMsg = []
        if self.outputInFloat16.get() and not self.saveSumPowerSpectra.get():
            errorMsg.append("'Save sum of power spectra?' must be set to Yes when writing in float16.")
        return errorMsg

    # --------------------------- UTILS functions -----------------------------
    def getInFilesPath(self):
        return self._getExtraPath('inFiles')

    def getMotionCorrSubtomosCmd(self):
        gainFile = self.inputTiltSeriesM.get().getGain()
        cmd = '--use_own --j 1 '
        cmd += f'--i {join(self.getInFilesPath(), IN_TS_STAR)} '
        cmd += f'--o {self._getExtraPath()} '
        if self.outputInFloat16.get():
            cmd += '--float16 '
        if self.saveEvenOdd.get():
            cmd += '--even_odd_split '
        if self.saveSumPowerSpectra.get():
            cmd += f'--grouping_for_ps {self.sumPowerSpectraNFrames.get()} '
        # MOTION TAB
        cmd += f'--bfactor {self.bFactor.get()} '
        cmd += f'--patch_x {self.patchX.get()} '
        cmd += f'--patch_y {self.patchY.get()} '
        cmd += f'--group_frames {self.groupFrames.get()} '
        cmd += f'--bin_factor {self.binningFactor.get()} '
        if gainFile:
            cmd += f'--gainref {gainFile} '
            cmd += f'--gain_rot {self.gainRot.get()} '
            cmd += f'--gain_flip {self.gainFlip.get()} '
        if self.defectFile.get():
            cmd += f'--defect_file {self.defectFile.get()} '
            cmd += f'--eer_grouping {self.eerFractionation.get()}'
        # EXTRA PARAMS
        cmd += self._genExtraParamsCmd()
        return cmd

