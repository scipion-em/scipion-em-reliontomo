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
import numpy as np

from pwem.protocols import EMProtocol
from pyworkflow import BETA
from pyworkflow.protocol import PointerParam, IntParam, GE, LE, PathParam
from pyworkflow.utils import Message
from reliontomo.constants import BOX_SIZE_VALS, OUT_TOMOS_STAR
from reliontomo.utils import getFileFromDataPrepProt


class ProtRelionPerParticlePerTiltBase(EMProtocol):
    """Base protocol used for the getting the frame alignment and ctf-refinment"""

    _devStatus = BETA

    _boxSize4Est = None

    # -------------------------- DEFINE param functions -----------------------

    def _defineParams(self, form):
        form.addSection(label=Message.LABEL_INPUT)
        form.addParam('inputPrepareDataProt', PointerParam,
                      pointerClass='ProtRelionPrepareData',
                      label="Data preparation protocol",
                      important=True,
                      allowsNull=False)
        form.addParam('inPseudoSubtomos', PointerParam,
                      pointerClass='SetOfPseudoSubtomograms',
                      label="Input pseudo-subtomograms",
                      important=True,
                      allowsNull=False)
        form.addParam('inputTrajectory', PathParam,
                      label="Particle trajectories set (optional)",
                      allowsNull=True)
        form.addParam('half1map', PathParam,
                      label='Path to half map 1',
                      help="Provide one of the two reference half-reconstructions MRC map files.")
        form.addParam('half2map', PathParam,
                      label='Path to half map 2',
                      help="Provide one of the two reference half-reconstructions MRC map files.")
        form.addParam('inRefMask', PointerParam,
                      pointerClass='Mask',
                      label="Input reference mask")
        form.addParam('inputPostProcess', PathParam,
                      label="Reference FSC (postprocess.star)",
                      allowsNull=True,
                      help="Input STAR file from a relion_postprocess job.")

    @staticmethod
    def _insertBoxSizeForEstimationParam(form):
        form.addParam('boxSize', IntParam,
                      label='Box size for estimation (pix)',
                      default=128,
                      allowsNull=False,
                      validators=[GE(32), LE(512)],
                      help="Box size to be used for the estimation. Note that this can be larger than the box size "
                           "of the reference map. A sufficiently large box size allows more of the high-frequency "
                           "signal to be captured that has been delocalized by the CTF.")

    # -------------------------- INSERT steps functions -----------------------
    def _insertAllSteps(self):
        pass

    # -------------------------- UTILS functions -----------------------------
    def _initialize(self):
        self._findClosestAdmittedVal()

    def _findClosestAdmittedVal(self):
        validVals = np.array(BOX_SIZE_VALS)
        # Find index of minimum value
        ind = np.where(validVals == np.amin(validVals - self.boxSize.get()))[0].tolist()[0]
        self._boxSize4Est = BOX_SIZE_VALS[ind]

    def _genIOCommand(self):
        cmd = '--p %s ' % self.inPseudoSubtomos.get().getStarFile()
        cmd += '--t %s ' % getFileFromDataPrepProt(self, OUT_TOMOS_STAR)
        cmd += '--o %s ' % self._getExtraPath()
        if self.inputTrajectory.get():
            cmd += '--mot %s ' % self.inputTrajectory.get()
        cmd += '--ref1 %s ' % self.half1map1.get()
        cmd += '--ref2 %s ' % self.half1map2.get()
        cmd += '--mask %s ' % self.inRefMask.get().getFileName()
        if self.inputPostProcess.get():
            cmd += '--fsc %s ' % self.inputPostProcess.get()







