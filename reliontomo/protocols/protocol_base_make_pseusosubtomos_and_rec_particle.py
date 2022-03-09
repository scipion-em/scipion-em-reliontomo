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
from pwem.protocols import EMProtocol
from pyworkflow import BETA
from pyworkflow.protocol import PointerParam, LEVEL_ADVANCED, IntParam, FloatParam
from pyworkflow.utils import Message


class ProtRelionMakePseudoSubtomoAndRecParticleBase(EMProtocol):
    """Reconstruct particle and make pseudo-subtomograms base class"""

    _label = None
    _devStatus = BETA

    # -------------------------- DEFINE param functions -----------------------
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.inParticlesStar = None
        self.inTomosStar = None
        self.inPrecedents = None

    def _defineParams(self, form):
        form.addSection(label=Message.LABEL_INPUT)

        form.addParam('inOptSet', PointerParam,
                      pointerClass='RelionParticles',
                      label='Input particles')
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

        form.addParallelSection(threads=1, mpi=1)

    # -------------------------- INSERT steps functions -----------------------
    def _insertAllSteps(self):
        pass

    # # -------------------------- INFO functions -------------------------------
    # def _validate(self):
    #     errorMsg = []
    #     if self.inputParticles.get():
    #         if not isPseudoSubtomogram(self.inputParticles.get().getFirstElement()):
    #             errorMsg.append('Introduced subtomograms do not contain the required data to be considered '
    #                             'pseudosubtomograms. This set can be generated using the output of the protocol for '
    #                             'preparing the data for Relion 4.')

    # # --------------------------- UTILS functions -----------------------------
    def _genCommonCmd(self):
        inRelionParticles = self.inOptSet.get()
        cmd = ''
        cmd += '--t %s ' % inRelionParticles.getTomograms()
        cmd += '--p %s ' % inRelionParticles.getParticles()
        if inRelionParticles.getTrajectories():
            cmd += '--mot %s ' % inRelionParticles.getTrajectories()
        cmd += '--b %i ' % self.boxSize.get()
        cmd += '--crop %i ' % self.croppedBoxSize.get()
        cmd += '--bin %.1f ' % self.binningFactor.get()
        cmd += '--SNR %.2f ' % self.snrWiener.get()
        cmd += '--j %i ' % self.numberOfThreads.get()
        return cmd

    def getNewSamplingRate(self):
        """It will be the tilt series sampling rate (bin 1) multiplied by the binning factor introduced"""
        return self.inOptSet.get().getTsSamplingRate() * self.binningFactor.get()


