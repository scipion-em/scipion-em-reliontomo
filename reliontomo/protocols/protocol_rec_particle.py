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
from pyworkflow.protocol import StringParam
from reliontomo import Plugin
from reliontomo.protocols.protocol_make_pseusosubtomos_and_rec_particle_base import \
    ProtRelionMakePseudoSubtomoAndRecParticleBase
from reliontomo.utils import getProgram, genSymmetryTable


class ProtRelionReconstructParticle(ProtRelionMakePseudoSubtomoAndRecParticleBase):
    """Reconstruct particle"""

    _label = 'Reconstruct particle'

    def __init__(self, **args):
        ProtRelionMakePseudoSubtomoAndRecParticleBase.__init__(self, **args)

    # -------------------------- DEFINE param functions -----------------------

    def _defineParams(self, form):
        ProtRelionMakePseudoSubtomoAndRecParticleBase._defineParams(self, form)
        form.addSection(label='Reconstruct particle')
        form.addParam('symmetry', StringParam,
                      label='Symmetry group',
                      default='C1',
                      help='Symmetry libraries have been copied from XMIPP. As such, with the exception of tetrahedral '
                           'symmetry, they comply with '
                           'https://relion.readthedocs.io/en/latest/Reference/Bibliography.html#id23. '
                           'Possible values [notation label] are described below:\n\n'
                           '%s' % json.dumps(genSymmetryTable(), indent=1))

    # -------------------------- INSERT steps functions -----------------------
    def _insertAllSteps(self):
        self._insertFunctionStep(self._relionReconstructParticle)

    # -------------------------- STEPS functions ------------------------------
    def _relionReconstructParticle(self):
        Plugin.runRelionTomo(self, 'relion_tomo_reconstruct_particle', self._genRecParticleCmd(),
                             numberOfMpi=self.numberOfMpi.get())

    # # -------------------------- INFO functions -------------------------------
    def _validate(self):
        pass

    # # --------------------------- UTILS functions -----------------------------
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