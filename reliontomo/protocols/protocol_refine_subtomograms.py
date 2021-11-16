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
import glob
import re

from emtable import Table

from pwem import ALIGN_PROJ
from pwem.convert.headers import fixVolume
from pwem.objects import FSC
from pyworkflow import BETA
from pyworkflow.object import Integer
from reliontomo.convert import convert40_tomo
from reliontomo.objects import SetOfPseudoSubtomograms
from reliontomo.protocols.protocol_base_refine import ProtRelionRefineBase
from reliontomo import Plugin
from os import listdir
from os.path import isfile, join, getmtime
from pyworkflow.protocol import PointerParam, LEVEL_ADVANCED, FloatParam, StringParam, BooleanParam, EnumParam
from pyworkflow.utils import moveFile
from reliontomo.constants import ANGULAR_SAMPLING_LIST, SYMMETRY_HELP_MSG
from reliontomo.utils import getProgram
from tomo.objects import AverageSubTomogram
from tomo.protocols import ProtTomoBase


class ProtRelionRefineSubtomograms(ProtRelionRefineBase, ProtTomoBase):
    """Auto-refinement of subtomograms."""

    _label = 'Auto-refinement of subtomograms'
    _devStatus = BETA
    FILE_KEYS = ['data', 'optimiser', 'sampling']
    PREFIXES = ['half1_', 'half2_']

    def __init__(self, **args):
        ProtRelionRefineBase.__init__(self, **args)
        self.reader = None

    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        self._defineInputParams(form)
        self._defineReferenceParams(form)
        self._defineCTFParams(form)
        self._defineOptimisationParamsCommon2All(form)
        self._defineAutoSamplingParams(form)
        self._defineComputeParams(form)
        ProtRelionRefineBase._defineAdditionalParams(form)

    @staticmethod
    def _defineInputParams(form):
        ProtRelionRefineBase._defineIOParams(form)
        form.addParam('referenceVolume', PointerParam,
                      pointerClass='Volume',
                      allowsNull=False,
                      label='Reference volume',
                      help='Initial reference 3D map, it should have the same dimensions and the same '
                           'pixel size as your input particles.')
        form.addParam('solventMask', PointerParam,
                      pointerClass='VolumeMask',
                      label='Reference mask (optional)',
                      allowsNull=True,
                      help='A volume mask containing a (soft) mask with the same dimensions '
                           'as the reference(s), and values between 0 and 1, with 1 being 100% protein '
                           'and 0 being 100% solvent. The reconstructed reference map will be multiplied '
                           'by this mask. If no mask is given, a soft spherical mask based on the <radius> '
                           'of the mask for the experimental images will be applied.\n\n'
                           'In some cases, for example for non-empty icosahedral viruses, it is also useful '
                           'to use a second mask. Check _Advaced_ parameters to select another volume mask.')
        form.addParam('solventMask2', PointerParam,
                      pointerClass='VolumeMask',
                      expertLevel=LEVEL_ADVANCED,
                      allowsNull=True,
                      label='Second reference mask (optional)',
                      help='For all white (value 1) pixels in this second mask the '
                           'corresponding pixels in the reconstructed map are set to the average value of '
                           'these pixels. Thereby, for example, the higher density inside the virion may be '
                           'set to a constant. Note that this second mask should have one-values inside the '
                           'virion and zero-values in the capsid and the solvent areas.')

    @staticmethod
    def _defineReferenceParams(form):
        form.addSection(label='Reference')
        form.addParam('isMapAbsoluteGreyScale', BooleanParam,
                      default=True,
                      label='Is initial 3D map on absolute greyscale?',
                      help='Perform CC-calculation in the first iteration (use this if references are not on the '
                           'absolute intensity scale). See detailed explanation below:\n\n '
                           'Probabilities are calculated based on a Gaussian noise model,'
                           'which contains a squared difference term between the reference and the experimental '
                           'image.\n\n This has a consequence that the reference needs to be on the same absolute '
                           'intensity greyscale as the experimental images. RELION and XMIPP reconstruct maps at '
                           'their absolute intensity greyscale. Other packages may perform internal normalisations of '
                           'the reference density, which will result in incorrect grey-scales. But, if the map was'
                           'reconstructed in RELION or in XMIPP, set this option to Yes, otherwise set it to No.\n\n'
                           'If set to No, RELION will use a (grey-scale invariant) cross-correlation criterion in the '
                           'first iteration, and prior to the second iteration the map will be filtered again using '
                           'the initial low-pass filter. This procedure is relatively quick and typically does not '
                           'negatively affect the outcome of the subsequent map refinement. Therefore, if in doubt it '
                           'is recommended to set this option to No.')
        form.addParam('initialLowPassFilterA', FloatParam,
                      default=30,
                      label='Initial low-pass filter (A)',
                      help='It is recommended to strongly low-pass filter your initial reference map. '
                           'If it has not yet been low-pass filtered, it may be done internally using this option. '
                           'If set to 0, no low-pass filter will be applied to the initial reference(s).')
        ProtRelionRefineBase._insertSymmetryParam(form)

    @staticmethod
    def _defineOptimisationParamsCommon2All(form):
        ProtRelionRefineSubtomograms._insertOptimisationSection(form)
        ProtRelionRefineBase._insertMaskDiameterParam(form)
        ProtRelionRefineBase._insertZeroMaskParam(form)
        ProtRelionRefineBase._insertFlattenSolventParam(form)
        form.addParam('solventCorrectFSC', BooleanParam,
                      default=False,
                      condition='solventMask',
                      label='Correct FSC curve for the effects of the solvent mask?',
                      help="If set to Yes, then instead of using unmasked maps to calculate the gold-standard FSCs "
                           "during refinement, masked half-maps are used and a post-processing-like correction of "
                           "the FSC curves (with phase-randomisation) is performed every iteration.\n\n"
                           "This only works when a reference mask is provided. This may yield "
                           "higher-resolution maps, especially when the mask contains only a relatively small "
                           "volume inside the box.")

    @staticmethod
    def _defineAutoSamplingParams(form):
        form.addSection(label='Auto-sampling')
        ProtRelionRefineBase._insertAngularCommonParams(form)
        form.addParam('localSearchAutoSampling', EnumParam,
                      default=4,
                      choices=ANGULAR_SAMPLING_LIST,
                      label='Local searches from auto-sampling',
                      help="Minimum healpix order (before oversampling) from which autosampling procedure will "
                           "use local searches.\n\n"
                           "In the automated procedure to increase the angular samplings, local angular "
                           "searches of -6/+6 times the sampling rate will be used from this angular sampling rate "
                           "onwards. For most lower-symmetric particles a value of 1.8 degrees will be sufficient. "
                           "Perhaps icosahedral symmetries may benefit from a smaller value such as 0.9 degrees.")
        ProtRelionRefineSubtomograms._insertRelaxSymmetry(form)
        form.addParam('useFinerAngularSampling', BooleanParam,
                      default=False,
                      label='Use finer angular sampling faster?',
                      help="If set to Yes, then let auto-refinement proceed faster with finer angular samplings. "
                           "Two additional conditions will be considered:\n\n "
                           "\t-Angular sampling will go down despite changes still happening in the angles.\n"
                           "\t-Angular sampling will go down if the current resolution already requires that sampling\n"
                           "\t at the edge of the particle.\n\nThis option will make the computation faster, but "
                           "hasn't been tested for many cases for potential loss in reconstruction quality upon "
                           "convergence.")

    @staticmethod
    def _defineComputeParams(form):
        ProtRelionRefineBase._defineComputeParams(form)
        form.addParam('skipPadding', BooleanParam,
                      default=False,
                      label='Skip padding?',
                      help="If set to Yes, the calculations will not use padding in Fourier space for better "
                           "interpolation in the references. Otherwise, references are padded 2x before Fourier "
                           "transforms are calculated. Skipping padding (i.e. use --pad 1) gives nearly as good "
                           "results as using --pad 2, but some artifacts may appear in the corners from signal "
                           "that is folded back.")

    @staticmethod
    def _insertRelaxSymmetry(form, condition=True):
        form.addParam('relaxSym', StringParam,
                      allowsNull=True,
                      condition=condition,
                      label='Symmetry to be relaxed',
                      help="With this option, poses related to the standard local angular search range by the given "
                           "point group will also be explored. For example, if you have a pseudo-symmetric dimer A-A', "
                           "refinement or classification in C1 with symmetry relaxation by C2 might be able to improve "
                           "distinction between A and A'. Note that the reference must be more-or-less aligned to the "
                           "convention of (pseudo-)symmetry operators. For details, see Ilca et al 2019 and Abrishami "
                           "et al 2020 cited in the About dialog.\n\n%s" % SYMMETRY_HELP_MSG)

    # -------------------------- INSERT steps functions -----------------------
    def _insertAllSteps(self):
        self._initialize()
        self._insertFunctionStep(self._autoRefine)
        self._insertFunctionStep(self.createOutputStep)

    # -------------------------- STEPS functions ------------------------------
    def _initialize(self):
        """ This function is mean to be called after the
        working dir for the protocol have been set. (maybe after recovery from mapper)
        """
        self._createFilenameTemplates()
        self._createIterTemplates()

    def _createFilenameTemplates(self):
        """ Centralize how files are called for iterations and references. """
        self.extraIter = self._getExtraPath('_it%(iter)03d_')
        myDict = {
            'input_star': self._getPath('input_particles.star'),
            'input_mrcs': self._getPath('input_particles.mrcs'),
            'data_scipion': self.extraIter + 'data_scipion.sqlite',
            'projections': self.extraIter + '%(half)sclass%(ref3d)03d_projections.sqlite',
            'classes_scipion': self.extraIter + 'classes_scipion.sqlite',
            'data': self.extraIter + 'data.star',
            'model': self.extraIter + 'model.star',
            'optimiser': self.extraIter + 'optimiser.star',
            'angularDist_xmipp': self.extraIter + 'angularDist_xmipp.xmd',
            'all_avgPmax': self._getPath('iterations_avgPmax.star'),
            'all_changes': self._getPath('iterations_changes.star'),
            'selected_volumes': self._getPath('selected_volumes_xmipp.xmd'),
            'dataFinal': self._getExtraPath("_data.star"),
            'modelFinal': self._getExtraPath("_model.star"),
            'finalvolume': self._getExtraPath("_class%(ref3d)03d.mrc"),
            'final_half1_volume': self._getExtraPath("_half1_class%(ref3d)03d_unfil.mrc"),
            'final_half2_volume': self._getExtraPath("_half2_class%(ref3d)03d_unfil.mrc"),
            'finalSGDvolume': self._getExtraPath("r_it%(iter)03d_class%(ref3d)03d.mrc"),
            'preprocess_particles': self._getPath("preprocess_particles.mrcs"),
            'preprocess_particles_star': self._getPath("preprocess_particles.star"),
            'preprocess_particles_prefix': "preprocess_particles"
        }
        # add to keys, data.star, optimiser.star and sampling.star
        for key in self.FILE_KEYS:
            myDict[key] = self.extraIter + '%s.star' % key
            key_xmipp = key + '_xmipp'
            myDict[key_xmipp] = self.extraIter + '%s.xmd' % key
        # add other keys that depends on prefixes
        for p in self.PREFIXES:
            myDict['%smodel' % p] = self.extraIter + '%smodel.star' % p
            myDict['%svolume' % p] = self.extraIter + p + 'class%(ref3d)03d.mrc'

        self._updateFilenamesDict(myDict)

    def _createIterTemplates(self):
        """ Setup the regex on how to find iterations. """
        self._iterTemplate = self._getFileName('data', iter=0).replace('000', '???')
        # Iterations will be identify by _itXXX_ where XXX is the iteration number
        # and is restricted to only 3 digits.
        self._iterRegex = re.compile('_it(\d{3})_')

    def _autoRefine(self):
        nMpi = self.numberOfMpi.get()
        Plugin.runRelionTomo(self, getProgram('relion_refine', nMpi), self._genAutoRefineCommand(), numberOfMpi=nMpi)

    def createOutputStep(self):
        subtomoSet = self.inputPseudoSubtomos.get()
        vol = AverageSubTomogram()
        volName = self._getExtraPath('_class001.mrc')
        fixVolume(volName)  # Fix header for xmipp to consider it a volume instead of a stack
        vol.setFileName(volName)
        vol.setSamplingRate(subtomoSet.getSamplingRate())
        pattern = '*it*half%s_class*.mrc'
        half1 = self._getLastFileName(self._getExtraPath(pattern % 1))
        half2 = self._getLastFileName(self._getExtraPath(pattern % 2))
        vol.setHalfMaps([half1, half2])

        outSubtomoSet = self._createSet(SetOfPseudoSubtomograms, 'pseudosubtomograms%s.sqlite', '')
        outSubtomoSet.copyInfo(subtomoSet)
        self.reader = convert40_tomo.Reader(alignType=ALIGN_PROJ, pixelSize=outSubtomoSet.getSamplingRate())
        self._fillDataFromIter(outSubtomoSet, self._lastIter())

        self._defineOutputs(outputVolume=vol)
        self._defineSourceRelation(subtomoSet, vol)
        self._defineOutputs(outputParticles=outSubtomoSet)
        self._defineTransformRelation(subtomoSet, outSubtomoSet)

        fsc = FSC(objLabel=self.getRunName())
        fn = self._getExtraPath("_model.star")
        table = Table(fileName=fn, tableName='model_class_1')
        resolution_inv = table.getColumnValues('rlnResolution')
        frc = table.getColumnValues('rlnGoldStandardFsc')
        fsc.setData(resolution_inv, frc)

        self._defineOutputs(outputFSC=fsc)
        self._defineSourceRelation(vol, fsc)

    # -------------------------- INFO functions -------------------------------
    def _validate(self):
        pass

    # --------------------------- UTILS functions -----------------------------
    def _genAutoRefineCommand(self):
        cmd = self._genBaseCommand()
        cmd += ' --auto_refine --split_random_halves --low_resol_join_halves 40 --norm --scale '
        # I/O args
        cmd += '--ref %s ' % self.referenceVolume.get().getFileName()
        if self.solventMask.get():
            cmd += '--solvent_mask %s ' % self.solventMask.get().getFileName()
        if self.solventMask2.get():
            cmd += '--solvent_mask2 %s ' % self.solventMask2.get().getFileName()
        # Reference args
        if self.isMapAbsoluteGreyScale.get():
            cmd += '--firstiter_cc '
        if self.initialLowPassFilterA.get():
            cmd += '--ini_high %.2f ' % self.initialLowPassFilterA.get()
        cmd += '--sym %s ' % self.symmetry.get()
        # Optimisation args
        if self.solventCorrectFSC.get():
            cmd += '--solvent_correct_fsc '
        if self.zeroMask.get():
            cmd += '--zero_mask '
        if self.flattenSolvent.get():
            cmd += '--flatten_solvent '
        # Angular sampling args
        cmd += '--healpix_order %i ' % self.angularSamplingDeg.get()  # - self.oversampling.get())
        cmd += '--offset_range %i ' % self.offsetSearchRangePix.get()
        cmd += '--offset_step %i ' % (self.offsetSearchStepPix.get() ** self.oversampling.get())
        cmd += '--auto_local_healpix_order %i ' % self.localSearchAutoSampling.get()  # - self.oversampling.get())
        if self.relaxSym.get():
            cmd += '--relax_sym %s ' % self.relaxSym.get()
        if self.useFinerAngularSampling.get():
            cmd += '--auto_ignore_angles --auto_resol_angles '
        # Compute args
        cmd += '--pad %i ' % (1 if self.skipPadding.get() else 2)

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

    @staticmethod
    def _getLastFileName(pattern):
        files = glob.glob(pattern)
        files.sort(key=getmtime)
        return files[-1]

    def _fillDataFromIter(self, subtomoClassesSet, iteration):
        tableName = 'particles'
        dataStar = tableName + '@' + self._getFileName('data', iter=iteration)
        for item in subtomoClassesSet:
            item.setAlignmentProj()

        mdIter = Table.iterRows(dataStar, key='rlnTomoParticleId')
        subtomoClassesSet.copyItems(self.inputPseudoSubtomos.get(),
                                    doClone=False,
                                    updateItemCallback=self._updateParticle,
                                    itemDataIterator=mdIter)

    def _updateParticle(self, particle, row):
        self.reader.setParticleTransform(particle, row)
        if not hasattr(particle, '_randomSubset'):
            particle._randomSubset = Integer()
            particle._randomSubset.set(row.rlnRandomSubset)

    def _getIterNumber(self, index):
        """ Return the list of iteration files, give the iterTemplate. """
        result = -1
        files = sorted(glob.glob(self._iterTemplate))
        if files:
            f = files[index]
            s = self._iterRegex.search(f)
            if s:
                result = int(s.group(1))  # group 1 is 3 digits iteration number
        return result

    def _lastIter(self):
        return self._getIterNumber(-1)
