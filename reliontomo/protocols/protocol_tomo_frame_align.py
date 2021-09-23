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
from pyworkflow.protocol import IntParam,  BooleanParam, GE, LE
from reliontomo.protocols.protocol_base_per_part_per_tilt import ProtRelionPerParticlePerTiltBase
from tomo.protocols import ProtTomoBase


class ProtRelionTomoFrameAlign(ProtRelionPerParticlePerTiltBase, ProtTomoBase):
    """Tomo frame align"""

    _label = 'Tomo frame align'

    def __init__(self, **args):
        ProtRelionPerParticlePerTiltBase.__init__(self, **args)

    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        ProtRelionPerParticlePerTiltBase._defineParams(self, form)
        form.addSection(label='Polish')
        ProtRelionPerParticlePerTiltBase._insertBoxSizeForEstimationParam(form)
        form.addParam('maxPosErr', IntParam,
                      label='Max position error (pix)',
                      default=5,
                      allowsNull=False,
                      validators=[GE(0), LE(64)],
                      help="maximal assumed error in the initial 2D particle-positions (distances between the "
                           "projected 3D positions and their true positions in the images), given in pixels.")
        form.addParam('doFlexAlign', BooleanParam,
                      label="Allow flexible alignment?",
                      default=False,
                      help="If set to No, only an optimal rigid shift will be applied to each frame (no iterative "
                           "optimisation).")
        form.addParam('doGlobalRigidAlign', BooleanParam,
                      label="Do global rigid shift alignment?",
                      default=False,
                      condition='not doFlexAlign',
                      help="If set to Yes, it estimates the rigid shift by aligning only the particles instead of by "
                           "predicting the entire micrographs.")

    # -------------------------- INSERT steps functions -----------------------
    def _insertAllSteps(self):
        pass
        # self._insertFunctionStep(self._relionMakePseudoSubtomos)
        # self._insertFunctionStep(self.createOutputStep)

    # -------------------------- STEPS functions ------------------------------

    # # -------------------------- INFO functions -------------------------------
    def _validate(self):
        pass

    # # --------------------------- UTILS functions -----------------------------


