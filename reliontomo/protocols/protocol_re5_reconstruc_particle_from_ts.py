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
from enum import Enum
from os.path import exists

from pwem.convert.headers import fixVolume
from pyworkflow.protocol import StringParam, FloatParam
from reliontomo import Plugin
from reliontomo.constants import SYMMETRY_HELP_MSG, IN_PARTICLES_STAR
from reliontomo.protocols.protocol_re5_base_extract_subtomos_and_rec_particle import \
    ProtRelion5ExtractSubtomoAndRecParticleBase
from tomo.objects import AverageSubTomogram


class outputObjects(Enum):
    average = AverageSubTomogram
    # postProcessVolume = VolumeMask
    # relionParticles = RelionSetOfPseudoSubtomograms


class ProtRelion5ReconstructParticle(ProtRelion5ExtractSubtomoAndRecParticleBase):
    """Reconstructs/averages from the tilt series projected particles"""

    _label = 'Reconstruct particle re5'
    _possibleOutputs = outputObjects

    def __init__(self, **args):
        super().__init__(**args)

    # -------------------------- DEFINE param functions -----------------------

    def _defineParams(self, form):
        super()._defineParams(form)
        form.addSection(label='Average')
        super()._defineCommonRecParams(form)
        form.addParam('symmetry', StringParam,
                      label='Symmetry group',
                      default='C1',
                      help=SYMMETRY_HELP_MSG)
        # form.addParam('solventMask', PointerParam,
        #               pointerClass='VolumeMask',
        #               label='FSC solvent mask (opt.)',
        #               allowsNull=True,
        #               help='Provide a soft mask to automatically estimate the postprocess FSC.')
        form.addParam('snrWiener', FloatParam,
                      label='Apply a Wiener filter with this SNR',
                      default=0,
                      help='If set to a positive value, apply a Wiener filter with this signal-to-noise ratio. If '
                           'omitted, the reconstruction will use a heuristic to prevent divisions by excessively '
                           'small numbers. Please note that using a low (even though realistic) SNR might wash out the '
                           'higher frequencies, which could make the map unsuitable to be used for further refinement.')
        self._defineExtraParams(form)

    # -------------------------- INSERT steps functions -----------------------
    def _insertAllSteps(self):
        self._insertFunctionStep(self.convertInputStep)
        self._insertFunctionStep(self.relionReconstructParticle)
        # if self.solventMask.get():
        #     self._insertFunctionStep(self.relionTomoMaskReference)
        self._insertFunctionStep(self.createOutputStep)

    # -------------------------- STEPS functions ------------------------------
    def relionReconstructParticle(self):
        cmd = self._genRecParticleCmd()
        try:
            Plugin.runRelionTomo(self, 'relion_tomo_reconstruct_particle_mpi', cmd,
                                 numberOfMpi=self.numberOfMpi.get())
        except:
            # The --mem argument should also be set using around 80-90% to keep a safety margin
            Plugin.runRelionTomo(self, 'relion_tomo_reconstruct_particle_mpi', cmd + '--mem 50 ',
                                 numberOfMpi=self.numberOfMpi.get())

    # def relionTomoMaskReference(self):
    #     Plugin.runRelionTomo(self, 'relion_tomo_make_reference', self._genTomoMaskRefCmd(),
    #                          numberOfMpi=self.numberOfMpi.get())

    def createOutputStep(self):
        inParticles = self.inReParticles.get()
        currentSamplingRate = inParticles.getTsSamplingRate() * self.binningFactor.get()
        halves = [self._getExtraPath('half1.mrc'), self._getExtraPath('half2.mrc')]

        # Fix headers to be interpreted as volumes instead of stacks
        [fixVolume(mrcFile) for mrcFile in glob.glob(self._getExtraPath('*.mrc'))]

        # Output average
        vol = AverageSubTomogram()
        vol.setFileName(self._getExtraPath('merged.mrc'))
        if exists(halves[0]):
            vol.setHalfMaps(halves)
        vol.setSamplingRate(currentSamplingRate)
        self._defineOutputs(**{outputObjects.average.name: vol})

        # # Output solvent mask
        # if self.solvenself._getExtraPath(IN_TOMOS_STAR)}tMask.get():
        #     postProccesMrc = self._genPostProcessOutputMrcFile(POST_PROCESS_MRC)
        #     postProccesMrc.setHalfMaps(halves)
        #     postProccesMrc.setSamplingRate(currentSamplingRate)
        #     self._defineOutputs(**{outputObjects.postProcessVolume.name: postProccesMrc})
        #     self._defineSourceRelation(inParticles, postProccesMrc)

        self._defineSourceRelation(inParticles, vol)

        # # Create the output set with the new optimization set
        # outParticles = self.genRelionParticles(boxSize=self.boxSize.get(),
        #                                        binningFactor=self.binningFactor.get())
        # self._defineOutputs(**{outputObjects.average.name: outParticles})

    # -------------------------- INFO functions -------------------------------
    def _validate(self):
        validateMsg = super()._validate()
        if not validateMsg:
            if self.numberOfMpi.get() == 1:
                validateMsg.append('The number of MPI must be greater than 1')
        return validateMsg

    # --------------------------- UTILS functions -----------------------------
    def _genRecParticleCmd(self):
        cmd = [
            f'--p {self._getExtraPath(IN_PARTICLES_STAR)}',
            f'--t {self.inReParticles.get().getTomogramsStar()}',
            f'--o {self._getExtraPath()}',
            self._genCommonExtractAndRecCmd(),
            f'--sym {self.symmetry.get()}',
            # Note:
            #   --j: number of threads used for the non-reconstruction parts of the program (e.g. symmetry application
            #        or gridding correction). This should be set to the number of CPU cores available.
            #   --j_out: number of threads that compute partial reconstructions in parallel. This is faster, but it
            #        requires additional memory for each thread. When used together with the --mem argument, this number
            #        will be reduced to (approximately) maintain the imposed memory limitation.
            #   --j_in: number of threads to be used for each partial reconstruction. This is a slower way to
            #        parallelise the procedure, but it does not require additional memory. Unless memory is limited,
            #        the --j_out option should be preferred. The product of --j_out and --j_in should not exceed the
            #        number of CPU cores available.
            f'--j_out {self.numberOfThreads.get()}',
            '--j_in 1'
        ]
        if self.snrWiener.get() > 0:
            cmd.append(f'--SNR {self.snrWiener.get():.2f}')
        return ' '.join(cmd)

    # def _genTomoMaskRefCmd(self):
    #     inParticles = self.inReParticles.get()
    #     cmd = ''
    #     cmd += '--t %s ' % inParticles.getTomogramsStar()
    #     cmd += '--p %s ' % self.getOutStarFileName()
    #     cmd += '--rec %s ' % self._getExtraPath()
    #     cmd += '--o %s ' % self._getExtraPath()
    #     # cmd += '--mask %s ' % self.solventMask.get().getFileName()
    #     cmd += '--angpix %.2f ' % (inParticles.getTsSamplingRate() * self.binningFactor.get())
    #     return cmd
