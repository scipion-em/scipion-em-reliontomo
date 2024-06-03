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

from pyworkflow.protocol import PointerParam, IntParam, FloatParam, BooleanParam, EnumParam, StringParam
from pyworkflow.utils import Message
from reliontomo import Plugin
from reliontomo.constants import IN_TS_STAR, TOMOGRAMS_DIR
from reliontomo.convert import convert50_tomo
from reliontomo.protocols.protocol_base_relion import ProtRelionTomoBase
from tomo.objects import SetOfTomograms

# Reconstruct options
SINGLE_TOMO = 0
ALL_TOMOS = 1


class tomoRecOtputObjects(Enum):
    tomograms = SetOfTomograms


class ProtRelion5TomoReconstruct(ProtRelionTomoBase):
    """Tomograms reconstruction"""

    _label = 'Tomograms reconstruction'
    _possibleOutputs = tomoRecOtputObjects

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

        # -------------------------- DEFINE param functions -----------------------

    def _defineParams(self, form):
        form.addSection(label=Message.LABEL_INPUT)
        form.addParam('inTsSet', PointerParam,
                      pointerClass='SetOfTiltSeries',
                      important=True,
                      label='Tilt-series')
        form.addParam('unbinnedWidth', IntParam,
                      default=-1,
                      allowsNull=False,
                      important=True,
                      label='Unbinned tomogram width (px)',
                      help="The tomogram X-dimension in unbinned pixels. It is recommended to use "
                           "a slightly larger tomogram volume than the actual size of the images so that if "
                           "they rotate, all pixels will still be in the tomogram.\n"
                           "If set to -1, the corresponding image dimension with an extra 10% will be used, e. g., "
                           "and image with X-dim equal to 4000 px would be passed to Relion as width = 4400 px.")
        form.addParam('unbinnedHeight', IntParam,
                      default=-1,
                      allowsNull=False,
                      important=True,
                      label='Unbinned tomogram height (px)',
                      help="The tomogram Y-dimension in unbinned pixels. See the help of parameter 'Unbinned tomogram "
                           "width (px)' for more details.")
        form.addParam('unbinnedThickness', IntParam,
                      default=-1,
                      allowsNull=False,
                      important=True,
                      label='Unbinned tomogram thickness (px)',
                      help="The tomogram Z-dimension in unbinned pixels. For your own data, you may want to test a "
                           "few values here to make sure the tomogram thickness is not too small to contain your "
                           "entire sample. If you intend to denoise your tomograms later with cryoCARE, it is better "
                           "not to pick a tomogram thickness that is much greater than the thickness of your sample, "
                           "because the denoising protocol randomly extracts subtomograms from your tomograms and you "
                           "don’t want too many without signal.\n"
                           "If set to -1, a thickness of 5/16 of the image width will be passed to Relion.")
        form.addParam('binnedPixSize', FloatParam,
                      allowsNull=False,
                      imortant=True,
                      label='Binned pixel size (Å/px)',
                      help='The tomogram will be downscaled to this pixel size. Typically, the larger the pixel size, '
                           'the faster the tomogram reconstruction and the less space the tomograms occupy on disk.')
        form.addParam('recTomoMode', EnumParam,
                      display=EnumParam.DISPLAY_HLIST,
                      choices=['Single tomogram', 'All tomograms'],
                      default=ALL_TOMOS,
                      label='Choose a reconstruction option. If the option Single tomograms is selected, then the '
                            'program will only reconstruct the tomogram chosen by tilt-series id.')
        form.addParam('tomoId', StringParam,
                      condition='recTomoMode == %s' % SINGLE_TOMO,
                      label='Tomogram to be reconstructed')
        form.addParam('genEvenOddTomos', BooleanParam,
                      default=False,
                      label='Generate the odd/even tomograms?',
                      help="Generate the odd/even tomograms that can be used for example to denoise with cryoCARE. For "
                           "this option to work, the introduced tilt-series must have their even/odd tilt-series "
                           "generated in the motion correction step.")
        form.addParam('tiltAngleOffset', FloatParam,
                      default=0,
                      label='Tilt angle offset (deg',
                      help="The tomogram tilt angles will all be changed by this value. This may be useful to "
                           "reconstruct lamellae that are all milled under a given angle. All tomograms will be "
                           "reconstructed with the same offset.")
        # TODO: add the params related to the 2D sums of the central Z-slices?

    # -------------------------- INSERT steps functions -----------------------
    def _insertAllSteps(self):
        self.convertInputStep()

    # -------------------------- STEPS functions ------------------------------
    def convertInputStep(self):
        outPath = self._getExtraPath()
        writer = convert50_tomo.Writer()
        # Aligned tilt-series star files: the one corresponding to the set and each TS star file
        writer.alignedTsSet2Star(self.inTsSet.get(), outPath)

    def reconstructTomogramsStep(self):
        Plugin.runRelionTomo(self, 'relion_tomo_reconstruct_tomogram', self.genTomoRecCmd())

    # --------------------------- UTILS functions -----------------------------
    def genTomoRecCmd(self):
        cmd = [
            f'--t {self._getExtraPath(IN_TS_STAR)}',
            f'--o {self._getExtraPath(TOMOGRAMS_DIR)}',
            f'--w {self.unbinnedWidth.get()}',
            f'--h {self.unbinnedHeight.get()}',
            f'--d {self.unbinnedThickness.get()}',
            f'--binned_angpix {self.binnedPixSize.get:.3f}',
            f'--j {self.numberOfThreads.get()}'
        ]
        if self.genEvenOddTomos.get():
            cmd.append('----generate_split_tomograms')
        if self.recTomoMode.get() == SINGLE_TOMO:
            cmd.append(f'--tn {self.tomoId.get()}')
        return ' '.join(cmd)

    # -------------------------- INFO functions -------------------------------
    def _validate(self):
        errorMsg = []
        if self.genEvenOddTomos.get() and not self.inTsSet.get().hasOddEven():
            errorMsg.append('Odd/even tomograms cannot be generated as the introduced tilt-series do not have odd/even '
                            'tilt-series associated.')
        return errorMsg

    def _warnings(self):
        warnMsg = []
        if not self.inTsSet.get().hasAlignment():
            warnMsg.append('The introduced tilt-series seems not to have alignment data.')
        return warnMsg

    def _summary(self):
        summary = []
        if self.isFinished():
            if self.recTomoMode.get() == SINGLE_TOMO:
                summary.append('The selected tomogram was *%s*.' % self.tomoId.get())
        return summary
