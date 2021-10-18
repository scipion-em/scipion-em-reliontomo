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
from os import mkdir, listdir
from os.path import join, exists
from imod.utils import generateDefocusIMODFileFromObject
from pwem.protocols import EMProtocol
from pyworkflow import BETA
from pyworkflow.protocol import PointerParam, PathParam, BooleanParam, LEVEL_ADVANCED, EnumParam
from reliontomo import Plugin
from reliontomo.constants import IN_TOMOS_STAR, IN_SUBTOMOS_STAR, OUT_TOMOS_STAR
from reliontomo.convert import writeSetOfTomograms, writeSetOfSubtomograms

# eTomo data source choices
ETOMO_FROM_PROT = 0
ETOMO_FROM_DIR = 1

# Other constants
DEFOCUS = 'defocus'


class ProtRelionPrepareDataJJ(EMProtocol):
    """Prepare data for Relion 4
    """
    _label = 'JJ Prepare data for Relion 4'
    _devStatus = BETA
    tomoSet = None
    acquisition = None
    inTomosStar = 'inTomos.star'
    inSubtomosStar = 'inSubtomos.star'

    # -------------------------- DEFINE param functions -----------------------

    def _defineParams(self, form):
        form.addSection(label='CTFTomoSeries')

        form.addParam('inputCtfTs', PointerParam,
                      pointerClass='SetOfCTFTomoSeries',
                      label="Input CTF tomo series",
                      important=True,
                      allowsNull=False,
                      help='Select the input set of CTF tomo series from the project.')

    # -------------------------- INSERT steps functions -----------------------
    def _insertAllSteps(self):
        self._insertFunctionStep(self._convertInputStep)

    # -------------------------- STEPS functions ------------------------------

    def _convertInputStep(self):
        # Generate defocus files
        for ctfTomo in self.inputCtfTs.get():
            defocusPath = self._getExtraPath(DEFOCUS, ctfTomo.getTsId())
            mkdir(self._getExtraPath(DEFOCUS))
            mkdir(defocusPath)
            generateDefocusIMODFileFromObject(ctfTomo,
                                              join(defocusPath, ctfTomo.getTsId() + '.' + DEFOCUS),
                                              isRelion=True)
