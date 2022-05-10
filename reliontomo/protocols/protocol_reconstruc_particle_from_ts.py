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
from pyworkflow.protocol import StringParam, PointerParam
from reliontomo import Plugin
from reliontomo.constants import SYMMETRY_HELP_MSG, OPTIMISATION_SET_STAR
from reliontomo.objects import relionTomoMetadata
from reliontomo.protocols.protocol_base_make_pseusosubtomos_and_rec_particle import \
    ProtRelionMakePseudoSubtomoAndRecParticleBase
from tomo.objects import AverageSubTomogram


class outputObjects(Enum):
    relionParticles = relionTomoMetadata
    average = AverageSubTomogram


class ProtRelionReconstructParticle(ProtRelionMakePseudoSubtomoAndRecParticleBase):
    """Reconstruct particle from the original tilt series images"""

    _label = 'Reconstruct particle from tilt series'
    _possibleOutputs = outputObjects

    def __init__(self, **args):
        ProtRelionMakePseudoSubtomoAndRecParticleBase.__init__(self, **args)

    # -------------------------- DEFINE param functions -----------------------

    def _defineParams(self, form):
        ProtRelionMakePseudoSubtomoAndRecParticleBase._defineParams(self, form)
        form.addSection(label='Reconstruct particle')
        form.addParam('symmetry', StringParam,
                      label='Symmetry group',
                      default='C1',
                      help=SYMMETRY_HELP_MSG)
        form.addParam('solventMask', PointerParam,
                      pointerClass='VolumeMask',
                      label='FSC solvent mask (opt.)',
                      allowsNull=True,
                      help='Provide a soft mask to automatically estimate the postprocess FSC.')

    # -------------------------- INSERT steps functions -----------------------
    def _insertAllSteps(self):
        self._insertFunctionStep(self.relionReconstructParticle)
        if self.solventMask.get():
            self._insertFunctionStep(self.relionTomoMaskReference)
        self._insertFunctionStep(self.createOutputStep)

    # -------------------------- STEPS functions ------------------------------
    def relionReconstructParticle(self):
        cmd = self._genRecParticleCmd()
        try:
            Plugin.runRelionTomo(self, 'relion_tomo_reconstruct_particle', cmd, numberOfMpi=self.numberOfMpi.get())
        except:
            # The --mem argument should also be set using around 80-90% to keep a safety margin
            Plugin.runRelionTomo(self, 'relion_tomo_reconstruct_particle', cmd + '--mem 50 ',
                                 numberOfMpi=self.numberOfMpi.get())

    def relionTomoMaskReference(self):
        Plugin.runRelionTomo(self, 'relion_tomo_make_reference', self._genTomoMaskRefCmd(),
                             numberOfMpi=self.numberOfMpi.get())

    def createOutputStep(self):
        inOptSet = self.inOptSet.get()
        # Output Relion Particles
        relionParticles = relionTomoMetadata(optimSetStar=self._getExtraPath(OPTIMISATION_SET_STAR),
                                             tsSamplingRate=inOptSet.getTsSamplingRate(),
                                             relionBinning=self.binningFactor.get(),
                                             nParticles=inOptSet.getNumParticles())

        # Output average
        vol = AverageSubTomogram()
        vol.setFileName(self._getExtraPath('merged.mrc'))
        vol.setHalfMaps([self._getExtraPath('half1.mrc'), self._getExtraPath('half2.mrc')])
        vol.setSamplingRate(relionParticles.getCurrentSamplingRate())

        self._defineOutputs(**{outputObjects.average.name: vol,
                               outputObjects.relionParticles.name: relionParticles})
        self._defineSourceRelation(inOptSet, relionParticles)
        self._defineSourceRelation(inOptSet, vol)

    # -------------------------- INFO functions -------------------------------
    def _validate(self):
        pass

    # --------------------------- UTILS functions -----------------------------
    def _genRecParticleCmd(self):
        cmd = self._genCommonCmd()
        cmd += '--o %s ' % self._getExtraPath()
        cmd += '--sym %s ' % self.symmetry.get()
        # Note:
        #   --j: number of threads used for the non-reconstruction parts of the program (e.g. symmetry application
        #        or gridding correction). This should be set to the number of CPU cores available.
        #   --j_out: number of threads that compute partial reconstructions in parallel. This is faster, but it
        #        requires additional memory for each thread. When used together with the --mem argument, this number
        #        will be reduced to (approximately) maintain the imposed memory limitation.
        #   --j_in: number of threads to be used for each partial reconstruction. This is a slower way to parallelise
        #        the procedure, but it does not require additional memory. Unless memory is limited, the --j_out option
        #        should be preferred. The product of --j_out and --j_in should not exceed the number of CPU cores
        #        available.
        cmd += '--j_out %i ' % self.numberOfThreads.get()
        cmd += '--j_in %i ' % 1
        return cmd

    def _genTomoMaskRefCmd(self):
        optSet = self.inOptSet.get()
        cmd = ''
        cmd += '--t %s ' % optSet.getTomograms()
        cmd += '--p %s ' % optSet.getParticles()
        cmd += '--rec %s ' % self._getExtraPath()
        cmd += '--o %s ' % self._getExtraPath()
        cmd += '--mask %s ' % self.solventMask.get().getFileName()

        return cmd
