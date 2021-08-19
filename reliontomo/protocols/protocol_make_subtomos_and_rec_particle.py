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
from pwem.protocols import EMProtocol
from pyworkflow import BETA
from pyworkflow.protocol import PointerParam, LEVEL_ADVANCED, IntParam, FloatParam, StringParam, BooleanParam
from pyworkflow.utils import Message
from reliontomo import Plugin
from reliontomo.constants import OUT_TOMOS_STAR, OUT_SUBTOMOS_STAR, REC_PARTICLES_DIR, PSEUDO_SUBTOMOS_DIR


class ProtRelionMakePseudoSubtomoAndRecParticle(EMProtocol):
    """Generate initial particle model"""
    _label = 'Generate initial particle model'
    _devStatus = BETA

    # -------------------------- DEFINE param functions -----------------------

    def _defineParams(self, form):
        form.addSection(label=Message.LABEL_INPUT)

        form.addParam('inputPrepareDataProt', PointerParam,
                      pointerClass='ProtRelionPrepareData',
                      label="Data preparation protocol",
                      important=True,
                      allowsNull=False)
        group = form.addGroup('Reconstruct')
        group.addParam('boxSize', IntParam,
                       label='Box size (pix.)',
                       important=True,
                       allowsNull=False,
                       help='Box size, in pixels,  of the reconstruction. Note that this is independent of the '
                            'box size used to refine the particle. This allows the user to construct a 3D map of '
                            'arbitrary size to gain an overview of the structure surrounding the particle. A '
                            'sufficiently large box size also allows more of the high-frequency signal to be '
                            'captured that has been delocalized by the CTF.')
        group.addParam('croppedBoxSize', IntParam,
                       label="Cropped box size (pix.)",
                       allowsNull=True,
                       help='Cropped box size in pixels. If set, the program will output an additional set of '
                            'maps that have been cropped to this size. This is useful if a map is desired that '
                            'is smaller than the box size required to retrieve the CTF-delocalized signal.')
        group.addParam('binningFactor', FloatParam,
                       label='Binning factor',
                       default=1,
                       allowsNull=False,
                       help='Downsampling (binning) factor. Note that this does not alter the box size. The '
                            'reconstructed region instead becomes larger.')
        group.addParam('snrWiener', FloatParam,
                       label='Apply a Wiener filter with this SNR',
                       default=-1,
                       expertLevel=LEVEL_ADVANCED,
                       help='Assumed signal-to-noise ratio (negative means use a heuristic to prevent divisions by '
                            'excessively small numbers.) Please note that using a low (even though realistic) SNR '
                            'might wash out the higher frequencies, which could make the map unsuitable to be used '
                            'for further refinement. More information about the Wiener Filter can be found here: '
                            'https://en.wikipedia.org/wiki/Wiener_filter')

        form.addSection(label='Reconstruct pseudo-Subtomograms')
        form.addParam('applyConeWeight', BooleanParam,
                      label='Apply cone weight?',
                      default=False,
                      help='Down weight a cone in Fourier space along the Z axis (as defined by the coordinate system '
                           'of the particle). This is useful for particles embedded in a membrane, as it can prevent '
                           'the alignment from being driven by the membrane signal (the signal of a planar membrane is '
                           'localised within one line in 3D Fourier space). Note that the coordinate system of a '
                           'particle is given by both the subtomogram orientation (if defined) and the particle '
                           'orientation. This allows the user to first obtain a membrane-driven alignment, and to then '
                           'specifically suppress the signal in that direction.')
        form.addParam('coneAngle', FloatParam,
                      label='Cone angle (deg.)',
                      condition='applyConeWeight',
                      default=10,
                      help='It is the (full) opening angle of the cone to be suppressed, given in degrees. This angle '
                           'should  include both the uncertainty about the membrane orientation and its variation '
                           'across the region represented in the subtomogram.')
        form.addParam('coneWidthZ0', IntParam,
                      label='Cone width at Z = 0 (pix.)',
                      condition='applyConeWeight',
                      default=2,
                      expertLevel=LEVEL_ADVANCED)
        form.addParam('keepParticlesDark', BooleanParam,
                      label='Do not invert contrast (keep particles dark)',
                      default=False,
                      expertLevel=LEVEL_ADVANCED)
        form.addParam('write3DCtfs', BooleanParam,
                      label='Write 3D CTFs?',
                      default=False,
                      expertLevel=LEVEL_ADVANCED)
        form.addParam('writeCtfCorrectedSubtomos', BooleanParam,
                      label='Write CTF-corrected subtomograms?',
                      default=False,
                      expertLevel=LEVEL_ADVANCED)

        form.addSection(label='Reconstruct particle')
        form.addParam('symmetry', StringParam,
                      label='Symmetry group',
                      default='C1',
                      help='Symmetry libraries have been copied from XMIPP. As such, with the exception of tetrahedral '
                           'symmetry, they comply with '
                           'https://relion.readthedocs.io/en/latest/Reference/Bibliography.html#id23. '
                           'Possible values [notation label] are described below:\n\n'
                           '%s' % json.dumps(self._genSymmetryTable(), indent=1))

        form.addParallelSection(threads=0, mpi=1)

    # -------------------------- INSERT steps functions -----------------------
    def _insertAllSteps(self):
        self._insertFunctionStep(self._relionReconstructParticle)
        self._insertFunctionStep(self._relionMakePseudoSubtomos)

    # -------------------------- STEPS functions ------------------------------
    def _relionReconstructParticle(self):
        self.runJob(self._getProgram('relion_tomo_reconstruct_particle'),
                    self._genRecParticleCmd(),
                    env=Plugin.getEnviron())

    def _relionMakePseudoSubtomos(self):
        self.runJob(self._getProgram('relion_tomo_subtomo'),
                    self._genMakePseudoSubtomoCmd(),
                    env=Plugin.getEnviron())

    # # -------------------------- INFO functions -------------------------------
    def _validate(self):
        pass

    # # --------------------------- UTILS functions -----------------------------
    def _getProgram(self, program):
        """ Get the program name depending on the MPI use or not. """
        if self.numberOfMpi.get() > 1:
            program += '_mpi'
        return program

    def _genCommonCmd(self):
        cmd = ''
        cmd += '--t %s ' % self._getFileFromDataPrepProt(OUT_TOMOS_STAR)
        cmd += '--p %s ' % self._getFileFromDataPrepProt(OUT_SUBTOMOS_STAR)
        cmd += '--b %i ' % self.boxSize.get()
        cmd += '--crop %i ' % self.croppedBoxSize.get()
        cmd += '--bin %.1f ' % self.binningFactor.get()
        cmd += '--SNR %.2f ' % self.snrWiener.get()
        cmd += '--j %i ' % self.numberOfMpi.get()
        return cmd

    def _genRecParticleCmd(self):
        cmd = self._genCommonCmd()
        cmd += '--o %s ' % self._getExtraPath(REC_PARTICLES_DIR)
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
        cmd += '--j_out %i ' % self.numberOfMpi.get()
        cmd += '--j_in %i ' % 1
        return cmd

    def _genMakePseudoSubtomoCmd(self):
        cmd = self._genCommonCmd()
        cmd += '--o %s ' % self._getExtraPath(PSEUDO_SUBTOMOS_DIR)
        if self.applyConeWeight.get():
            cmd += '--cone_weight '
            cmd += '--cone_angle %.2f ' % self.coneAngle.get()
            cmd += '--cone_sig0 %i ' % self.coneWidthZ0.get()
        if self.keepParticlesDark.get():
            cmd += ' --no_ic '
        if self.write3DCtfs.get():
            cmd += '--ctf '
        if self.writeCtfCorrectedSubtomos.get():
            cmd += '--div '
        return cmd

    def _getFileFromDataPrepProt(self, fileName):
        return self.inputPrepareDataProt.get()._getExtraPath(fileName)

    @staticmethod
    def _genSymmetryTable():
        jsonData = '[{"group": "Asymmetric", "notation": "C1", "origin": "User-defined", "orientation": "User-defined"},' \
                   '{"group": "Cyclic", "notation": "C<n>", "origin": "On symm axis, Z user-defined", "orientation": "Symm axis on Z"},' \
                   '{"group": "Dihedral", "notation": "D<n>", "origin": "Intersection of symm axes", "orientation": "principle symm axis on Z, 2-fold on X"},' \
                   '{"group": "Tetrahedral", "notation": "T", "origin": "Intersection of symm axes", "orientation": "3-fold axis on Z (deviating from Heymann et al!)"},' \
                   '{"group": "Octahedral", "notation": "O", "origin": "Intersection of symm axes", "orientation": "4-fold axes on X, Y, Z"},' \
                   '{"group": "Icosahedral", "notation": "I<n>", "origin": "Intersection of symm axes", "orientation": "**"}]'

        return json.loads(jsonData)
