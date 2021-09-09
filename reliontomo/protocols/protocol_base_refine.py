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
from os import remove
from os.path import abspath, exists
from pwem.protocols import EMProtocol
from pyworkflow import BETA
from pyworkflow.protocol import PointerParam, LEVEL_ADVANCED, IntParam, StringParam, BooleanParam, \
    EnumParam, PathParam
from pyworkflow.utils import Message
from reliontomo.constants import ANGULAR_SAMPLING_LIST, OUT_SUBTOMOS_STAR, SYMMETRY_HELP_MSG


class ProtRelionRefineBase(EMProtocol):
    """Base protocol used for the getting the initial model and performing the auto-refinment"""

    _devStatus = BETA

    # -------------------------- DEFINE param functions -----------------------

    def _defineParams(self, form):
        form.addParallelSection(threads=1, mpi=1)

    @staticmethod
    def _defineIOParams(form):
        form.addSection(label=Message.LABEL_INPUT)
        form.addParam('inputPseudoSubtomosProt', PointerParam,
                      pointerClass='ProtRelionMakePseudoSubtomograms',
                      label="Data preparation protocol",
                      important=True,
                      allowsNull=False)

    @staticmethod
    def _defineCTFParams(form):
        form.addSection(label='CTF')
        form.addParam('doCTF', BooleanParam,
                      default=True,
                      label='Do CTF-correction?',
                      help='If set to Yes, CTFs will be corrected inside the MAP refinement. '
                           'The resulting algorithm intrinsically implements the optimal linear, '
                           'or Wiener filter. Note that input particles should contains CTF parameters.')
        form.addParam('ignoreCTFUntilFirstPeak', BooleanParam,
                      default=False,
                      label='Ignore CTFs until first peak?',
                      help='If set to Yes, then CTF-amplitude correction will only be performed from the first peak '
                           'of each CTF onward. This can be useful if the CTF model is inadequate at the lowest '
                           'resolution. Still, in general using higher amplitude contrast on the CTFs (e.g. 10-20%) '
                           'often yields better results. Therefore, this option is not generally recommended.')
        form.addParam('ctfPhaseFlipped', BooleanParam,
                      default=False,
                      label='Has the data been CTF phase-flipped?',
                      expertLevel=LEVEL_ADVANCED)
        form.addParam('padCtf', BooleanParam,
                      default=False,
                      label='Perform CTF padding to treat CTF aliaising better?',
                      expertLevel=LEVEL_ADVANCED)
        form.addParam('ctfUncorrectedRef', BooleanParam,
                      default=False,
                      label='Have the input references not been CTF-amplitude corrected?',
                      expertLevel=LEVEL_ADVANCED)

    @staticmethod
    def _defineOptimisationParams(form):
        form.addSection(label='Optimisation')
        form.addParam('maskDiameter', IntParam,
                      allowsNull=False,
                      label='Circular mask diameter (Å)',
                      help='Diameter of the circular mask that will be applied to the experimental images '
                           '(in Angstroms)')
        form.addParam('zeroMask', BooleanParam,
                      allowsNull=False,
                      label='Mask surrounding background in particles to zero?',
                      default=False,
                      help='Diameter of the circular mask that will be applied to the experimental images '
                           '(in Angstroms)')

    @staticmethod
    def _defineComputeParams(form):
        form.addSection(label='Compute')
        form.addParam('noParallelDiscIO', BooleanParam,
                      default=False,
                      label='Do not let MPI processes access the disc simultaneously',
                      help='Do NOT let parallel (MPI) processes access the disc simultaneously (use '
                           'this option with NFS).')
        form.addParam('pooledSubtomos', IntParam,
                      default=1,
                      label='Number of pooled particles',
                      help='Number of images to pool for each thread task.')
        form.addParam('skipGridding', BooleanParam,
                      default=False,
                      label='Skip griding?',
                      help='Skip gridding in the Maximization step in the Expectation-Maximization algorithm. '
                           'If this option is set to Yes, more memory will be consumed during the protocol execution, '
                           'but it will be faster.')
        form.addParam('allParticlesRam', BooleanParam,
                      default=False,
                      label='Pre-read all particles into RAM?',
                      help='If set to Yes, the leader process read all particles into memory. Be careful you have '
                           'enough RAM for large data sets!')
        form.addParam('combineItersDisc', BooleanParam,
                      default=False,
                      label='Combine iterations through disc?',
                      help='If set to Yes, the large arrays of summed weights will be sent through the MPI network '
                           'instead of writing large files to disc.')
        form.addParam('scratchDir', PathParam,
                      label='Copy particles to scratch directory',
                      help='If provided, particle stacks will be copied to this local scratch disk prior for '
                           'refinement.')
        form.addParam('doGpu', BooleanParam,
                      default=False,
                      label='Use GPU acceleration?',
                      help='If set to Yes, it will use available gpu resources for some calculations.')
        form.addParam('gpusToUse', StringParam,
                      condition='doGpu',
                      default='0',
                      label='GPUs to use:',
                      help='It can be used to provide a list of which GPUs (e. g. "0:1:2:3") to use. MPI-processes are '
                           'separated by ":", threads by ",". For example: "0,0:1,1:0,0:1,1"')

    @staticmethod
    def _defineAdditionalParams(form):
        form.addSection(label='Additional')
        form.addParam('keepOnlyLastIterFiles',
                      BooleanParam,
                      default=False,
                      expertLevel=LEVEL_ADVANCED,
                      label="Keep only files from last iteration?",
                      help="If Yes is chosen, only the files which correspond to the last iteration will be saved "
                           "in the protocol's extra directory. Otherwise, files corresponding to each iteration "
                           "will be kept.")
        form.addParam('oversampling', IntParam,
                      default=1,
                      expertLevel=LEVEL_ADVANCED,
                      label="Over-sampling",
                      help="Adaptive oversampling order to speed-up calculations (0=no oversampling, 1=2x, 2=4x, etc)")
        form.addParam('extraParams', StringParam,
                      label='Additional arguments',
                      help="In this box command-line arguments may be provided that are not generated by the GUI. This "
                           "may be useful for testing developmental options and/or expert use of the program, e.g: \n"
                           "--verb 1\n"
                           "--pad 2\n")

    @staticmethod
    def _addSymmetryParam(form):
        form.addParam('symmetry', StringParam,
                      label='Symmetry group',
                      default='C1',
                      help=SYMMETRY_HELP_MSG)

    @staticmethod
    def _addAngularCommonParams(form):
        form.addParam('angularSamplingDeg', EnumParam,
                      default=2,
                      choices=ANGULAR_SAMPLING_LIST,
                      label='Angular sampling interval (deg)',
                      help='There are only a few discrete angular samplings possible because '
                           'we use the HealPix library to generate the sampling of the first '
                           'two Euler angles on the sphere. The samplings are approximate numbers '
                           'and vary slightly over the sphere.')
        form.addParam('offsetSearchRangePix', IntParam,
                      default=6,
                      label='Offset search range (pix.)',
                      help='Probabilities will be calculated only for translations in a circle '
                           'with this radius (in pixels). The center of this circle changes at '
                           'every iteration and is placed at the optimal translation for each '
                           'image in the previous iteration.')
        form.addParam('offsetSearchStepPix', IntParam,
                      default=2,
                      label='Offset search step (pix.)',
                      help='Translations will be sampled with this step-size (in pixels). '
                           'Translational sampling is also done using the adaptive approach. '
                           'Therefore, if adaptive=1, the translations will first be evaluated'
                           'on a 2x coarser grid.')

    # -------------------------- INSERT steps functions -----------------------
    def _insertAllSteps(self):
        pass

    # -------------------------- INFO functions -------------------------------
    def _validate(self):
        pass

    # --------------------------- UTILS functions -----------------------------
    def _genCommonCommand(self):
        cmd = ''
        cmd += '--i %s ' % self.inputPseudoSubtomosProt.get()._getExtraPath(OUT_SUBTOMOS_STAR)
        cmd += '--o %s ' % self._getExtraPath()  # If not, Relion will concatenate it directly converting the
        # last part of the path into a prefix of the resulting filename
        cmd += '--j %i ' % self.numberOfThreads
        # CTF args
        if self.doCTF.get():
            cmd += '--ctf '
        if self.ignoreCTFUntilFirstPeak.get():
            cmd += '--ctf_intact_first_peak '
        if self.ctfPhaseFlipped.get():
            cmd += '--ctf_phase_flipped '
        if self.padCtf.get():
            cmd += '--pad_ctf '
        if self.ctfUncorrectedRef.get():
            cmd += '--ctf_uncorrected_ref '

        # Optimisation args
        cmd += '--particle_diameter %i ' % self.maskDiameter.get()
        if self.zeroMask.get():
            cmd += '--zero_mask '

        # Compute args
        if self.noParallelDiscIO.get():
            cmd += '--no_parallel_disc_io '
        cmd += '--pool %i ' % self.pooledSubtomos.get()
        if self.allParticlesRam.get():
            cmd += '--preread_images '
        if self.skipGridding.get():
            cmd += '--skip_gridding '
        if not self.combineItersDisc.get():
            cmd += '--dont_combine_weights_via_disc '
        if self.scratchDir.get():
            cmd += '--scratch_dir %s ' % self.scratchDir.get()
        if self.doGpu.get():
            cmd += '--gpu %s ' % self.gpusToUse.get()

        # Additional args
        cmd += 'oversampling %i' % self.oversampling.get()
        if self.extraParams.get():
            cmd += ' ' + self.extraParams.get()

        return cmd

    def _applyKeepIterFilesUserSelection(self):
        if self.keepOnlyLastIterFiles.get():
            self._cleanUndesiredFiles()

    def _cleanUndesiredFiles(self):
        """Remove all files generated by relion_classify 3d excepting the ones which
        correspond to the last iteration. Example for iteration 25:
        relion_it025_class002.mrc
        relion_it025_class001.mrc
        relion_it025_model.star
        relion_it025_sampling.star
        relion_it025_optimiser.star
        relion_it025_data.star
        """
        itPref = 'relion_it'
        clPref = 'class'
        starExt = '.star'
        mrcExt = '.mrc'
        # Classify calculations related files
        calcFiles = ['data', 'model', 'optimiser', 'sampling']
        for i in range(self._lastIter()):
            for calcFile in calcFiles:
                fn = abspath(self._getExtraPath('{}{:03d}_{}{}'.format(
                    itPref, i, calcFile, starExt)))
                if exists(fn):
                    remove(fn)
            # Classes related files
            for itr in range(1, self.nClasses + 1):
                fn = abspath(self._getExtraPath('{}{:03d}_{}{:03d}{}'.format(
                    itPref, i, clPref, itr, mrcExt)))
                if exists(fn):
                    remove(fn)
