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
from reliontomo.protocols.protocol_base_refine import ProtRelionRefineBase
from tomo.objects import AverageSubTomogram
from reliontomo import Plugin
from os import listdir
from os.path import isfile, join
from pyworkflow.protocol import LEVEL_ADVANCED, IntParam, BooleanParam
from pyworkflow.utils import moveFile
from reliontomo.utils import getProgram


class ProtRelionDeNovoInitialModel(ProtRelionRefineBase):
    """Generate a de novo 3D initial model from the pseudo-subtomograms."""

    _label = 'Generate a de novo 3D initial model from the pseudo-subtomograms'

    def __init__(self, **args):
        ProtRelionRefineBase.__init__(self, **args)

    # -------------------------- DEFINE param functions -----------------------

    def _defineParams(self, form):
        ProtRelionRefineBase._defineParams(self, form)
        ProtRelionRefineBase._defineIOParams(form)
        ProtRelionRefineBase._defineCTFParams(form)
        self._defineOptimisationParams(form)
        ProtRelionRefineBase._defineComputeParams(form)
        ProtRelionRefineBase._defineAdditionalParams(form)

    @staticmethod
    def _defineOptimisationParams(form):
        ProtRelionRefineBase._defineOptimisationParams(form)
        form.addParam('maxNumberOfIterations', IntParam,
                      default=25,
                      label='Number of iterations',
                      help='Maximum number of iterations to be performed.')
        form.addParam('numberOfClasses', IntParam,
                      default=1,
                      label='Number of classes to be defined.')
        form.addParam('gradBasedOpt', BooleanParam,
                      default=False,
                      label='Perform gradient based optimisation',
                      expertLevel=LEVEL_ADVANCED,
                      help='Perform gradient based optimisation (instead of default expectation-maximization).')
        form.addParam('gradWriteIter', IntParam,
                      default=10,
                      label='Write out model every number of iterations',
                      expertLevel=LEVEL_ADVANCED,
                      help='Write out model every so many iterations during gradient refinement')
        form.addParam('noInitBlobs', BooleanParam,
                      default=False,
                      label='Switch off initializing models with random Gaussians?',
                      expertLevel=LEVEL_ADVANCED)
        form.addParam('flattenSolvent', BooleanParam,
                      default=False,
                      label='Flatten and enforce non-negative solvent?')
        ProtRelionRefineBase._addSymmetryParam(form)
        ProtRelionRefineBase._addAngularCommonParams(form)

    # -------------------------- INSERT steps functions -----------------------
    def _insertAllSteps(self):
        self._insertFunctionStep(self._generateDeNovo3DModel)
        self._insertFunctionStep(self.createOutputStep)

    # -------------------------- STEPS functions ------------------------------
    def _generateDeNovo3DModel(self):
        # Gradient based optimisation is not compatible with MPI (relion throws an exception mentioning it)
        nMpi = 1 if self.gradBasedOpt.get() else self.numberOfMpi.get()
        Plugin.runRelionTomo(self, getProgram('relion_refine', nMpi), self._genInitModelCommand(), numberOfMpi=nMpi)

    def createOutputStep(self):
        self._manageGeneratedFiles()
        vol = AverageSubTomogram()
        vol.setFileName(self._getExtraPath(self._getModelName()))
        vol.setSamplingRate(8.83)  # TODO: check how to get the sampling rate at this point of the pipeline
        self._defineOutputs(outputVolume=vol)

    # -------------------------- INFO functions -------------------------------
    def _validate(self):
        pass

    # --------------------------- UTILS functions -----------------------------
    def _genInitModelCommand(self):
        cmd = self._genCommonCommand()
        cmd += '--denovo_3dref '
        # Optimisation args
        cmd += '--iter %i ' % self.maxNumberOfIterations.get()
        cmd += '--K %i ' % self.numberOfClasses.get()
        if self.flattenSolvent.get():
            cmd += '--flatten_solvent '
        if self.gradBasedOpt.get():
            cmd += '--grad '
        if self.gradWriteIter.get():
            cmd += '--grad_write_iter %i ' % self.gradWriteIter.get()
        if self.noInitBlobs.get():
            cmd += '--no_init_blobs '
        cmd += '--sym %s ' % self.symmetry.get()
        cmd += '--healpix_order %i ' % self.angularSamplingDeg.get()
        cmd += '--offset_step %i ' % self.offsetSearchStepPix.get()
        cmd += '--offset_range %i ' % self.offsetSearchRangePix.get()
        return cmd

    def _getModelName(self):
        """generate the name of the volume following this pattern extra_it002_class001.mrc"""
        return 'it{:03d}_class001.mrc'.format(self.maxNumberOfIterations.get())

    def _manageGeneratedFiles(self):
        """There's some kind of bug in relion4 which makes it generate the file in the protocol base directory
        instead of the extra directory. It uses extra as a prefix of each generated file instead. Hence, until
        it's solved, the files will be moved to the extra directory and the prefix extra_ will be removed"""
        prefix = 'extra_'
        genFiles = [f for f in listdir(self._getPath()) if isfile(join(self._getPath(), f))]
        for f in genFiles:
            moveFile(self._getPath(f), self._getExtraPath(f.replace(prefix, '')))
