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
import json

from reliontomo.protocols.protocol_base_refine import ProtRelionRefineBase
from tomo.objects import AverageSubTomogram
from reliontomo import Plugin
from os import remove, listdir
from os.path import abspath, exists, isfile, join
from pwem.protocols import EMProtocol
from pyworkflow import BETA
from pyworkflow.protocol import PointerParam, LEVEL_ADVANCED, IntParam, FloatParam, StringParam, BooleanParam, \
    EnumParam, PathParam
from pyworkflow.utils import Message, moveFile
from reliontomo.constants import ANGULAR_SAMPLING_LIST, OUT_SUBTOMOS_STAR
from reliontomo.utils import genSymmetryTable, getProgram


class ProtRelionRefineSUbtomograms(ProtRelionRefineBase):
    """Subtomograms auto-refinement."""

    _label = 'Subtomograms auto-refinement'

    # -------------------------- INSERT steps functions -----------------------
    def _insertAllSteps(self):
        self._insertFunctionStep(self._autoRefine)
        self._insertFunctionStep(self.createOutputStep)

    # -------------------------- STEPS functions ------------------------------
    def _autoRefine(self):
        pass
        # # Gradient based optimisation is not compatible with MPI (relion throws an exception mentioning it)
        # nMpi = 1 if self.gradBasedOpt.get() else self.numberOfMpi.get()
        # Plugin.runRelionTomo(self, getProgram('relion_refine', nMpi), self._genCommand(), numberOfMpi=nMpi)

    def createOutputStep(self):
        pass

    # -------------------------- INFO functions -------------------------------
    def _validate(self):
        pass

    # --------------------------- UTILS functions -----------------------------
    def _genInitModelCommand(self):
        cmd = self._genCommand()
        cmd += '--denovo_3dref '
        return cmd

    def _getModelName(self):
        """generate the name of the volume following this pattern extra_it002_class001.mrc"""
        return 'it{:03d}_class001.mrc'.format(self.maxNumberOfIterations.get())

    def _manageGeneratedFiles(self):
        """There's some kind of bug in relion4 which makes it generate the file in the protocol base directory
        instead of the extra directory. It uses extra as a prefix of each generated file instead. Hence, until
        it's solved, the files will be moved to the extra directory and the prefix extra_ will be removed"""
        prefix = '_extra'
        genFiles = [f for f in listdir(self._getPath()) if isfile(join(self._getPath(), f))]
        for f in genFiles:
            moveFile(self._getPath(f), self._getExtraPath(f.replace(prefix, '')))
