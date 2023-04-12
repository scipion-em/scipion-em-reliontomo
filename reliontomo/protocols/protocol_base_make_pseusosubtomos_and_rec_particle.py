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
from pyworkflow import BETA
from pyworkflow.protocol import IntParam, FloatParam
from reliontomo.protocols.protocol_base_relion import ProtRelionTomoBase


class ProtRelionMakePseudoSubtomoAndRecParticleBase(ProtRelionTomoBase):
    """Reconstruct particle and make pseudo-subtomograms base class"""

    _label = None
    _devStatus = BETA

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        super()._defineCommonInputParams(form)
        form.addParallelSection(threads=1, mpi=1)

    @staticmethod
    def _defineCommonRecParams(form):
        form.addParam('boxSize', IntParam,
                      label='Box size (pix.)',
                      important=True,
                      allowsNull=False,
                      help='Box size, in pixels,  of the reconstruction. Note that this is independent of the '
                           'box size used to refine the particle. This allows the user to construct a 3D map of '
                           'arbitrary size to gain an overview of the structure surrounding the particle. A '
                           'sufficiently large box size also allows more of the high-frequency signal to be '
                           'captured that has been delocalized by the CTF.')
        form.addParam('croppedBoxSize', IntParam,
                      label="Cropped box size (pix.)",
                      allowsNull=True,
                      help='Cropped box size in pixels. If set, the program will output an additional set of '
                           'maps that have been cropped to this size. This is useful if a map is desired that '
                           'is smaller than the box size required to retrieve the CTF-delocalized signal.')
        form.addParam('binningFactor', FloatParam,
                      label='Binning factor',
                      default=1,
                      allowsNull=False,
                      help='Downsampling (binning) factor. Note that this does not alter the box size. The '
                           'reconstructed region instead becomes larger.')

    # -------------------------- INSERT steps functions -----------------------
    def _insertAllSteps(self):
        pass

    def convertInputStep(self):
        self.genInStarFile()

    def createOutputStep(self):
        return self.genRelionParticles(binningFactor=self.binningFactor.get(), boxSize=self.croppedBoxSize.get())

    # -------------------------- INFO functions -------------------------------

    # --------------------------- UTILS functions -----------------------------
    def _genCommonCmd(self):
        inRelionParticles = self.inReParticles.get()
        cmd = ''

        # Cancel this for now
        self.info("Using optimization_set: %s" % inRelionParticles.filesMaster)
        cmd += '--i %s ' % inRelionParticles.filesMaster

        # This would be either the particles start file in the set or a new generated one from the subtomo set.
        cmd += '--p %s ' % self.getOutStarFileName()

        if inRelionParticles.getTrajectories():
            cmd += '--mot %s ' % inRelionParticles.getTrajectories()


        cmd += '--b %i ' % self.boxSize.get()
        cmd += '--crop %i ' % self.croppedBoxSize.get()
        cmd += '--bin %.1f ' % self.binningFactor.get()
        cmd += '--j %i ' % self.numberOfThreads.get()
        return cmd

