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
from reliontomo.protocols import ProtRelionRefineSubtomograms, ProtRelionDeNovoInitialModel
from reliontomo.protocols.protocol_base_refine import ProtRelionRefineBase
from reliontomo import Plugin
from os import listdir
from os.path import isfile, join
from pyworkflow.protocol import PointerParam, LEVEL_ADVANCED, FloatParam, StringParam, BooleanParam, EnumParam, \
    LEVEL_NORMAL, GE, LE
from pyworkflow.utils import moveFile
from reliontomo.constants import ANGULAR_SAMPLING_LIST, SYMMETRY_HELP_MSG
from reliontomo.utils import getProgram


class ProtRelion3DClassifySubtomograms(ProtRelionRefineSubtomograms):
    """3D Classification of subtomograms."""

    _label = '3D Classification of subtomograms'

    def __init__(self, **args):
        ProtRelionRefineSubtomograms.__init__(self, **args)

    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        ProtRelionRefineSubtomograms._defineInputParams(form)
        ProtRelionRefineSubtomograms._defineReferenceParams(form)
        ProtRelionRefineSubtomograms._defineCTFParams(form)
        self._defineOptimisationParams(form)
        self._defineSamplingParams(form)
        ProtRelionRefineSubtomograms._defineComputeParams(form)
        ProtRelionRefineSubtomograms._defineAdditionalParams(form)

    @staticmethod
    def _defineOptimisationParams(form):
        ProtRelionRefineSubtomograms._insertOptimisationSection(form)
        ProtRelionRefineSubtomograms._insertNumOfClassesParam(form)
        ProtRelionRefineSubtomograms._insertRegularisationParam(form)
        ProtRelionRefineSubtomograms._insertNItersParam(form)
        form.addParam('useFastSubsets', BooleanParam,
                      label='Use fast subsets (for large data sets)?',
                      default=False,
                      help='If set to Yes, the first 5 iterations will be done with random subsets of only K*1500 '
                           'particles (K being the number of classes); the next 5 with K*4500 particles, the next '
                           '5 with 30% of the data set; and the final ones with all data. This was inspired by a '
                           'cisTEM implementation by Niko Grigorieff et al.')
        ProtRelionRefineSubtomograms._insertMaskDiameterParam(form)
        ProtRelionRefineSubtomograms._insertZeroMaskParam(form)
        form.addParam('limitResolutionEStepTo', FloatParam,
                      label='Limit resolution E-step to (Å)',
                      default=-1,
                      help='If set to a positive number, then the expectation step (i.e. the alignment) will be done '
                           'only including the Fourier components up to this resolution (in Angstroms).\n\nThis is '
                           'useful to prevent overfitting, as the classification runs in RELION are not guaranteed to '
                           'be 100% overfitting-free (unlike the 3D auto-refine with its gold-standard FSC). '
                           'In particular for very difficult data sets, e.g. of very small or featureless particles, '
                           'this has been shown to give much better class averages. \n\nIn such cases, values in the '
                           'range of 7-12 Angstroms have proven useful.')

    @staticmethod
    def _defineSamplingParams(form):
        form.addSection(label='Sampling')
        form.addParam('doImageAlignment', BooleanParam,
                      label='Perform image alignment?',
                      default=False,
                      help='If set to No, then rather than performing both alignment and classification, only '
                           'classification will be performed. This allows the use of very focused masks. It requires '
                           'that the optimal orientations of all particles are already stored in the input STAR file.')
        ProtRelionRefineBase._insertAngularCommonParams(form,
                                                        angSampling=2,
                                                        offsetRange=5,
                                                        offsetStep=1,
                                                        condition='doImageAlignment')
        form.addParam('doLocalAngleSearch', BooleanParam,
                      label='Perform local angular searches?',
                      default=False,
                      condition='doImageAlignment',
                      help="If set to Yes, then rather than performing exhaustive angular searches, local searches "
                           "within the range given below will be performed.\n\nA prior Gaussian distribution centered "
                           "at the optimal orientation in the previous iteration and with a stddev of 1/3 of the range "
                           "given below will be enforced.")
        form.addParam('localAngularSearchRange', FloatParam,
                      label='Local angular search range',
                      condition='doImageAlignment and doLocalAngleSearch',
                      default=5,
                      validators=[GE(0), LE(15)],
                      help="Local angular searches will be performed within +/- the given amount (in degrees) from "
                           "the optimal orientation in the previous iteration.\n\nA Gaussian prior (also see previous "
                           "option) will be applied, so that orientations closer to the optimal orientation in the "
                           "previous iteration will get higher weights than those further away.")
        ProtRelionRefineSubtomograms._insertRelaxSymmetry(form, condition='doImageAlignment and doLocalAngleSearch')

    # -------------------------- INSERT steps functions -----------------------
    def _insertAllSteps(self):
        self._insertFunctionStep(self._autoRefine)
        self._insertFunctionStep(self.createOutputStep)

    # -------------------------- STEPS functions ------------------------------
    def _autoRefine(self):
        nMpi = self.numberOfMpi.get()
        Plugin.runRelionTomo(self, getProgram('relion_refine', nMpi), self._genAutoRefineCommand(), numberOfMpi=nMpi)

    def createOutputStep(self):
        pass

    # -------------------------- INFO functions -------------------------------
    def _validate(self):
        pass

    # --------------------------- UTILS functions -----------------------------
    def _genAutoRefineCommand(self):
        cmd = self._genCommonCommand()
        cmd += '--auto_refine --split_random_halves --low_resol_join_halves 40 --norm --scale '
        # I/O args
        cmd += '--ref %s ' % self.referenceVolume.get()
        if self.solventMask.get():
            cmd += '--solvent_mask %s ' % self.solventMask.get()
        if self.solventMask2.get():
            cmd += '--solvent_mask2 %s ' % self.solventMask2.get()
        # Reference args
        if self.isMapAbsoluteGreyScale.get():
            cmd += 'firstiter_cc '
        if self.initialLowPassFilterA.get():
            cmd += '--ini_high %.2f ' % self.initialLowPassFilterA.get()
        # Optimisation args
        if self.solventCorrectFSC.get():
            cmd += '--solvent_correct_fsc '
        # Angular sampling args
        if self.localSearchAutoSampling.get():
            cmd += '--auto_local_healpix_order %i ' % self.localSearchAutoSampling.get()
        if self.relaxSym.get():
            cmd += '--relax_sym %s ' % self.relaxSym.get()
        if self.useFinerAngularSampling.get():
            cmd += '--auto_ignore_angles --auto_resol_angles '
        # Compute args
        cmd += '--pad %i' % (1 if self.skipPadding.get() else 2)

        return cmd

    def _getModelName(self):
        '''generate the name of the volume following this pattern extra_it002_class001.mrc'''
        return 'it{:03d}_class001.mrc'.format(self.maxNumberOfIterations.get())

    def _manageGeneratedFiles(self):
        '''There's some kind of bug in relion4 which makes it generate the file in the protocol base directory
        instead of the extra directory. It uses extra as a prefix of each generated file instead. Hence, until
        it's solved, the files will be moved to the extra directory and the prefix extra_ will be removed'''
        prefix = '_extra'
        genFiles = [f for f in listdir(self._getPath()) if isfile(join(self._getPath(), f))]
        for f in genFiles:
            moveFile(self._getPath(f), self._getExtraPath(f.replace(prefix, '')))
