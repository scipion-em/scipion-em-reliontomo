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
from os.path import exists
from typing import Union
from pwem.convert.headers import fixVolume
from pwem.objects import Volume, SetOfVolumes
from pyworkflow.object import Pointer
from reliontomo.constants import SYMMETRY_HELP_MSG
from reliontomo.protocols.protocol_base_refine import ProtRelionRefineBase
from reliontomo.protocols.protocol_base_relion import IS_RELION_50
from reliontomo import Plugin
from pyworkflow.protocol import LEVEL_ADVANCED
from reliontomo.utils import getProgram


class outputObjects(Enum):
    average = Volume
    averages = SetOfVolumes


class ProtRelionDeNovoInitialModel(ProtRelionRefineBase):
    """ Generate a de novo 3D initial model from the pseudo-subtomograms.

    This de novo 3D initial model allows to obtain a map without any prior knowledge.
    Provided you have a reasonable distribution of angular directions, this algorithm
    is likely to yield a suitable, low-resolution model that can subsequently be used
    for 3D classification or 3D auto-refine.\n

    Relion-4.0 uses a gradient-driven algorithm to generate a de novo 3D initial model
    from the pseudo-subtomograms.
    """

    _label = '3D initial model'
    _possibleOutputs = outputObjects

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    # -------------------------- DEFINE param functions -----------------------

    def _defineParams(self, form):
        self._defineIOParams(form)
        self._defineCTFParams(form)
        self._defineOptimisationParams(form)
        self._defineComputeParams(form)
        self._insertGpuParams(form)
        self._defineAdditionalParams(form)
        form.addParallelSection(threads=0, mpi=1)

    def _defineOptimisationParams(self, form):
        self._insertOptimisationSection(form)
        self._insertVdamMiniBatchesParam(form)
        self._insertRegularisationParam(form)
        self._insertNumOfClassesParam(form)
        self._insertMaskDiameterParam(form)
        self._insertFlattenSolventParam(form)
        helpDeNovo = 'The initial model is always generated in C1 and then aligned to and symmetrized ' \
                     'with the specified point group. If the automatic alignment fails, please manually  ' \
                     'rotate run_itNNN_class001.mrc (NNN is the number of iterations) so that it conforms ' \
                     'the symmetry convention.' + SYMMETRY_HELP_MSG
        self._insertSymmetryParam(form, helpDeNovo)
        self._insertDoInC1AndApplySymLaterParam(form)
        if IS_RELION_50:
            self._insertPriorWidthParam(form)
        self._insertAngularCommonParams(form,
                                        expertLevel=LEVEL_ADVANCED,
                                        angSampling=1,
                                        offsetRange=6,
                                        offsetStep=2)

    # -------------------------- INSERT steps functions -----------------------
    def _insertAllSteps(self):
        self._insertFunctionStep(self.convertInputStep, needsGPU=False)
        self._insertFunctionStep(self.generateDeNovo3DModel, needsGPU=True)
        self._insertFunctionStep(self.alignSymmetry, needsGPU=False)
        self._insertFunctionStep(self.createOutputStep, needsGPU=False)

    # -------------------------- STEPS functions ------------------------------
    def generateDeNovo3DModel(self):
        nMpi = self.numberOfMpi.get()
        Plugin.runRelionTomo(self,
                             getProgram('relion_refine', nMpi),
                             self._genInitModelCommand(),
                             numberOfMpi=nMpi)

    def alignSymmetry(self):
        if self.numberOfClasses.get() == 1:
            Plugin.runRelionTomo(self,
                                 'relion_align_symmetry',
                                 self._genApplySymCmd())
        else:
            for i in range(self.numberOfClasses.get()):
                Plugin.runRelionTomo(self,
                                     'relion_align_symmetry',
                                     self._genApplySymCmd(classIndex=i + 1))

    def createOutputStep(self):
        inRelionParticlesPointer = self.getInputParticles(returnPointer=True)
        if self.numberOfClasses.get() == 1:
            vol = self._createOutputModel(inRelionParticlesPointer)
            self._defineOutputs(**{outputObjects.average.name: vol})
            self._defineSourceRelation(inRelionParticlesPointer, vol)
        else:
            avgSet = SetOfVolumes.create(self._getPath(), template='averages%s.sqlite')
            avgSet.copyInfo(self.getInputParticles())
            for i in range(self.numberOfClasses.get()):
                vol = self._createOutputModel(inRelionParticlesPointer, classIndex=i + 1)
                avgSet.append(vol)

            self._defineOutputs(**{outputObjects.averages.name: avgSet})
            self._defineSourceRelation(inRelionParticlesPointer, avgSet)

    # -------------------------- INFO functions -------------------------------
    def _validate(self):
        errorMsg = []
        nMpi = self.numberOfMpi.get()
        if nMpi > 1:
            errorMsg.append('The initial volume can only run using 1 MPI.')

        if self.inReParticles.get().getBoxSize() % 2 != 0:
            errorMsg.append('The dimensions of the extracted pseudo-subtomograms '
                            'must be even!')
        return errorMsg

    # --------------------------- UTILS functions -----------------------------
    def _genInitModelCommand(self) -> str:
        # Common parameters from base protocol
        cmd = self._genBaseCommand(useOptimizationSet=False)

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
        if IS_RELION_50:
            cmd += '--sigma_tilt %i ' % self.priorWidthTiltAngle.get()
        return cmd

    def _genApplySymCmd(self, classIndex: Union[int, None] = None) -> str:
        cmd = '--o %s ' % self._getInitialModelOutFn(classIndex=classIndex)
        classIndex = 1 if classIndex is None else classIndex
        cmd += '--i %s ' % self._getExtraPath(self._getModelName(classIndex))
        if self.doInC1AndApplySymLater.get() and 'c1' not in self.symmetry.get().lower():
            cmd += '--sym %s ' % self.symmetry.get()
        else:
            cmd += '--sym C1 '
        cmd += '--apply_sym --select_largest_class '
        return cmd

    def _getModelName(self, classIndex: Union[int, None]) -> str:
        """generate the name of the volume following this pattern _it002_model.star"""
        classIndexStr = f'{classIndex:03d}' if classIndex is not None else ''
        return f'_it{self.nVdamMiniBatches.get():03d}_class{classIndexStr}.mrc'

    def _getInitialModelOutFn(self, classIndex: Union[int, None] = None) -> str:
        classIndexStr = f'_{classIndex:03d}' if classIndex is not None else ''
        return self._getExtraPath(f'initial_model{classIndexStr}.mrc')

    def _createOutputModel(self,
                           inRelionParticlesPointer: Pointer,
                           classIndex: Union[int, None] = None) -> Volume:
        vol = Volume()
        iniModelFile = self._getInitialModelOutFn(classIndex=classIndex)
        fixVolume(iniModelFile)  # Fix header to make it interpreted as volume instead of a stack by xmipp
        vol.setFileName(iniModelFile)
        vol.setSamplingRate(inRelionParticlesPointer.get().getSamplingRate())
        return vol
