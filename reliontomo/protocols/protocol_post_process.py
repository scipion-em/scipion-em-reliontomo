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
from os.path import join, dirname

import reliontomo
from pwem.protocols import EMProtocol
from pyworkflow import BETA
from pyworkflow.protocol import PointerParam, BooleanParam, FloatParam, EnumParam, GE, LE, IntParam
from pyworkflow.utils import Message
from reliontomo import Plugin
from reliontomo.constants import OUT_PARTICLES_STAR, COORD_X, COORD_Y, COORD_Z, SHIFTX_ANGST, SHIFTY_ANGST, \
    SHIFTZ_ANGST, ROT, TILT, PSI
from reliontomo.objects import relionTomoMetadata
from reliontomo.utils import genRelionParticles


resourcesPath = join(dirname(reliontomo.__file__), 'resources')


class outputObjects(Enum):
    outputRelionParticles = relionTomoMetadata


class detectoMftFiles(Enum):
    no_mft = ('No MFT file', None)
    de_20_300 = ('DE-20 at 300 KV', join(resourcesPath, 'DE-20_300KV.star'))
    falcon2_300 = ('Falcon II at 300 KV', join(resourcesPath, 'Falcon-II_300KV.star'))
    k2_summit_300 = ('K2-summit at 300 Kv', join(resourcesPath, 'K2-summit_300KV.star'))


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
        form.addParam('detectorMtf', EnumParam,
                      choices=[detectoMftFile.value(0) for detectoMftFile in detectoMftFiles],
                      default=detectoMftFiles.no_mft.value(0),
                      label='MTF of the detector',
                      help='If you know the MTF of your detector, provide it here. Curves provided were downloaded '
                           'from the RELION Wiki. Also see there for the exact format. If you do not know the MTF of '
                           'your detector and do not want to measure it, then by leaving this entry empty, you '
                           'include the MTF of your detector in your overall estimated B-factor upon sharpening the '
                           'map.Although that is probably slightly less accurate, the overall quality of your map will '
                           'probably not suffer very much.')
        form.addParam('origDetectorPixSize', FloatParam,
                      default=1,
                      validators=[GE(0.1), LE(2)],
                      condition='detectorMtf != %s' % detectoMftFiles.no_mft.value(0),
                      label='Original detector pixel size ((Å)/pix)',
                      help='This is the original pixel size (in Angstroms) in the raw (non-super-resolution!) '
                           'micrographs.')

    # -------------------------- INSERT steps functions -----------------------
    def _insertAllSteps(self):
        pass


    # def createOutputStep(self):
    #     # Output RelionParticles
    #     relionParticles = genRelionParticles(self._getExtraPath(), self.inOptSet.get())
    #     self._defineOutputs(**{outputObjects.outputRelionParticles.name: relionParticles})

    # -------------------------- INFO functions -------------------------------
    def _validate(self):
        pass

    # --------------------------- UTILS functions -----------------------------
