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
from pwem.convert.headers import fixVolume
from reliontomo.constants import INITIAL_MODEL
from reliontomo.protocols.protocol_base_refine import ProtRelionRefineBase
from tomo.objects import AverageSubTomogram
from reliontomo import Plugin
from pyworkflow.protocol import LEVEL_ADVANCED
from reliontomo.utils import getProgram


class outputObjects(Enum):
    average = AverageSubTomogram()


class ProtRelionDeNovoInitialModel(ProtRelionRefineBase):
    """Generate a de novo 3D initial model from the pseudo-subtomograms."""

    _label = 'De novo 3D initial model'
    _possibleOutputs = outputObjects

    def __init__(self, **args):
        ProtRelionRefineBase.__init__(self, **args)

    # -------------------------- DEFINE param functions -----------------------

    def _defineParams(self, form):
        ProtRelionRefineBase._defineIOParams(form)
        ProtRelionRefineBase._defineCTFParams(form)
        self._defineOptimisationParams(form)
        ProtRelionRefineBase._defineComputeParams(form)
        ProtRelionRefineBase._insertGpuParams(form)
        ProtRelionRefineBase._defineAdditionalParams(form)

    @staticmethod
    def _defineOptimisationParams(form):
        ProtRelionRefineBase._insertOptimisationSection(form)
        ProtRelionRefineBase._insertVdamMiniBatchesParam(form)
        ProtRelionRefineBase._insertRegularisationParam(form)
        ProtRelionRefineBase._insertNumOfClassesParam(form)
        ProtRelionRefineBase._insertMaskDiameterParam(form)
        ProtRelionRefineBase._insertFlattenSolventParam(form)
        ProtRelionRefineBase._insertSymmetryParam(form)
        ProtRelionRefineBase._insertDoInC1AndApplySymLaterParam(form)
        ProtRelionRefineBase._insertAngularCommonParams(form,
                                                        expertLevel=LEVEL_ADVANCED,
                                                        angSampling=1,
                                                        offsetRange=6,
                                                        offsetStep=2)

    # -------------------------- INSERT steps functions -----------------------
    def _insertAllSteps(self):
        self._insertFunctionStep(self.generateDeNovo3DModel)
        self._insertFunctionStep(self.alignSymmetry)
        self._insertFunctionStep(self.createOutputStep)

    # -------------------------- STEPS functions ------------------------------
    def generateDeNovo3DModel(self):
        Plugin.runRelionTomo(self, getProgram('relion_refine', self.numberOfMpi.get()), self._genInitModelCommand(),
                             numberOfMpi=self.numberOfMpi.get())

    def alignSymmetry(self):
        Plugin.runRelionTomo(self, 'relion_align_symmetry', self._genApplySymCmd())

    def createOutputStep(self):
        vol = AverageSubTomogram()
        iniModel = self._getExtraPath(INITIAL_MODEL)
        fixVolume(iniModel)  # Fix header to make it interpreted as volume instead of a stack by xmipp
        vol.setFileName(iniModel)
        vol.setSamplingRate(self.inOptSet.get().getCurrentSamplingRate())
        self._defineOutputs(**{outputObjects.average.name: vol})
        self._defineSourceRelation(self.inOptSet.get(), vol)

    # -------------------------- INFO functions -------------------------------

    # --------------------------- UTILS functions -----------------------------
    def _genInitModelCommand(self):
        # Common parameters from base protocol
        cmd = self._genBaseCommand()

        # Initial model specific commands
        cmd += ' --denovo_3dref --grad --zero_mask --auto_sampling --pad 1 '
        #   Optimisation args
        cmd += '--iter %i ' % self.nVdamMiniBatches.get()
        cmd += '--tau2_fudge %d ' % self.regularisation.get()
        cmd += '--K %i ' % self.numberOfClasses.get()
        if self.flattenSolvent.get():
            cmd += '--flatten_solvent '
        if self.doInC1AndApplySymLater.get():
            cmd += '--sym C1 '
        else:
            cmd += '--sym %s ' % self.symmetry.get()
        cmd += '--healpix_order %i ' % self.angularSamplingDeg.get()
        cmd += '--offset_step %i ' % self.offsetSearchStepPix.get()
        cmd += '--offset_range %d ' % self.offsetSearchRangePix.get()

        return cmd

    def _genApplySymCmd(self):
        cmd = '--i %s ' % self._getExtraPath(self._getModelName())
        cmd += '--o %s ' % self._getExtraPath(INITIAL_MODEL)
        if self.doInC1AndApplySymLater.get() and 'c1' not in self.symmetry.get().lower():
            cmd += '--sym %s ' % self.symmetry.get()
        else:
            cmd += '--sym C1 '
        cmd += '--apply_sym --select_largest_class '
        return cmd

    def _getModelName(self):
        """generate the name of the volume following this pattern _it002_model.star"""
        return '_it{:03d}_model.star'.format(self.nVdamMiniBatches.get())
