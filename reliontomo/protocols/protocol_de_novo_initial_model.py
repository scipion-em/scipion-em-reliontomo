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
from os import remove
from os.path import abspath, exists
from pwem.protocols import EMProtocol
from pyworkflow import BETA
from pyworkflow.protocol import PointerParam, LEVEL_ADVANCED, IntParam, FloatParam, StringParam, BooleanParam, \
    EnumParam, PathParam
from pyworkflow.utils import Message
from reliontomo.constants import ANGULAR_SAMPLING_LIST
from reliontomo.utils import genSymmetryTable


class ProtRelionDeNovoInitialModel(EMProtocol):
    """Generate a de novo 3D initial model from the pseudo-subtomograms."""

    _label = 'Generate a de novo 3D initial model from the pseudo-subtomograms'
    _devStatus = BETA

    # -------------------------- DEFINE param functions -----------------------

    def _defineParams(self, form):
        form.addSection(label=Message.LABEL_INPUT)
        form.addParam('inputPrepareDataProt', PointerParam,
                      pointerClass='ProtRelionPrepareData',
                      label="Data preparation protocol",
                      important=True,
                      allowsNull=False)

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

        form.addSection(label='Optimisation')
        form.addParam('numberOfIterations', IntParam,
                      default=25,
                      label='Number of iterations',
                      help='Number of iterations to be performed.')
        form.addParam('numberOfClasses', IntParam,
                      default=1,
                      label='Number of classes to be defined.')
        form.addParam('maskDiameter', FloatParam,
                      default=-1,
                      label='Circular mask diameter (Å)',
                      help='Diameter of the circular mask that will be applied to the experimental images '
                           '(in Angstroms)')
        form.addParam('flattenSolvent', BooleanParam,
                      default=False,
                      label='Flatten and enforce non-negative solvent?')
        form.addParam('symmetry', StringParam,
                      label='Symmetry group',
                      default='C1',
                      help='Symmetry libraries have been copied from XMIPP. As such, with the exception of tetrahedral '
                           'symmetry, they comply with '
                           'https://relion.readthedocs.io/en/latest/Reference/Bibliography.html#id23. '
                           'Possible values [notation label] are described below:\n\n'
                           '%s' % json.dumps(genSymmetryTable(), indent=1))
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

        form.addSection(label='Compute')
        form.addParam('noParallelDiscIO', BooleanParam,
                      default=False,
                      label='Use parallel disc I/O?',
                      help='Do NOT let parallel (MPI) processes access the disc simultaneously (use '
                           'this option with NFS).')
        form.addParam('pooledSubtomos', IntParam,
                      default=1,
                      label='Number of pooled particles',
                      help='Number of images to pool for each thread task.')
        form.addParam('allParticlesRam', BooleanParam,
                      default=False,
                      label='Pre-read all particles into RAM?',
                      help='If set to Yes, all particle images will be read into computer memory, which will greatly '
                           'speed up calculations on systems with slow disk access. However, one should of course be '
                           'careful with the amount of RAM available. Because particles are read in float-precision, '
                           'it will take \n( N * (box_size)^2 * 4 / (1024 * 1024 * 1024) ) Giga-bytes to read N '
                           'particles into RAM. For 100 thousand 200x200 images, that becomes 15Gb, or 60 Gb for the '
                           'same number of 400x400 particles. Remember that running a single MPI slave on each node '
                           'that runs as many threads as available cores will have access to all available RAM.\n\n'
                           'If parallel disc I/O is set to No, then only the master reads all particles into RAM and '
                           'sends those particles through the network to the MPI slaves during the refinement '
                           'iterations.')
        form.addParam('combineItersDisc', BooleanParam,
                      default=False,
                      label='Combine iterations through disc?',
                      help='If set to Yes, at the end of every iteration all MPI slaves will write out a large file '
                           'with their accumulated results. The MPI master will read in all these files, combine '
                           'them all, and write out a new file with the combined results. All MPI slaves will then '
                           'read in the combined results. This reduces heavy load on the network, but increases load '
                           'on the disc I/O. This will affect the time it takes between the progress-bar in the '
                           'expectation step reaching its end (the mouse gets to the cheese) and the start of the '
                           'ensuing maximisation step. It will depend on your system setup which is most efficient.')
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
                      label='GPUs to use:',
                      help='It can be used to provide a list of which GPUs (e. g. "0:1:2:3") to use. MPI-processes are '
                           'separated by ":", threads by ",". For example: "0,0:1,1:0,0:1,1"')

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

        form.addParallelSection(threads=0, mpi=1)

    # -------------------------- INSERT steps functions -----------------------
    def _insertAllSteps(self):
        pass

    # -------------------------- STEPS functions ------------------------------

    # # -------------------------- INFO functions -------------------------------
    def _validate(self):
        pass

    # # --------------------------- UTILS functions -----------------------------
    # if self.keepOnlyLastIterFiles:
    #     self._cleanUndesiredFiles()

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
