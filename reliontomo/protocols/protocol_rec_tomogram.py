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
from pwem.convert.headers import fixVolume
from pwem.protocols import EMProtocol
from pyworkflow import BETA
from pyworkflow.protocol.params import FloatParam, IntParam, StringParam, PointerParam, EnumParam
from pyworkflow.utils import Message
from reliontomo import Plugin
from reliontomo.protocols.protocol_prepare_data import outputObjects as prepareProtOutputs
from tomo.objects import Tomogram, SetOfTomograms
from tomo.utils import getObjFromRelation

# Reconstruct options
SINGLE_TOMO = 0
ALL_TOMOS = 1


class outputObjects(Enum):
    tomograms = SetOfTomograms


class ProtRelionTomoReconstruct(EMProtocol):
    """ This protocol reconstructs a single tomogram using Relion. It is very useful
    to check if the protocol "Prepare data" has been applied correctly (in terms of flip
    options, for example).
    """
    _label = 'Reconstruct tomograms from prepare data prot'
    _devStatus = BETA

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.inReParticles = None
        self.tomoSet = None
        self.tomoList = None
        self.outTomoSet = None
        self.outputSamplingRate = None

    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        form.addSection(label=Message.LABEL_INPUT)
        form.addParam('protPrepare', PointerParam,
                      pointerClass='ProtRelionPrepareData',
                      label='Prepare Data protocol',
                      help='It is very usefulto check if the protocol "Prepare data" has been applied correctly (in '
                           'terms of flip options, for example)')
        form.addParam('recTomoMode', EnumParam,
                      display=EnumParam.DISPLAY_HLIST,
                      choices=['Single tomogram', 'All tomograms'],
                      default=SINGLE_TOMO,
                      label='Choose a reconstruction option')
        form.addParam('tomoId', StringParam,
                      condition='recTomoMode == %s' % SINGLE_TOMO,
                      label='Tomogram to be reconstructed')
        form.addParam('binFactor', FloatParam,
                      label='Binning factor',
                      default=8,
                      help='The binning will be applied concerning the size of '
                           'the tomograms used for the picking.')
        group = form.addGroup('Tomogram shape (pix.)')
        group.addParam('width', IntParam,
                       default=-1,
                       label='Width',
                       help='If -1, the width considered will be of the original tilt series after having applied the '
                            'introduced binning factor.')
        group.addParam('height', IntParam,
                       default=-1,
                       label='Height',
                       help='If -1, the height considered will be of the original tilt series after having applied the '
                            'introduced binning factor.')
        group.addParam('thickness', IntParam,
                       default=-1,
                       label='Thickness',
                       help='If -1, the thickness considered will be of the original tilt series after having applied '
                            'the introduced binning factor.')

        form.addParallelSection(threads=4, mpi=0)

    # -------------------------- INSERT steps functions -----------------------
    def _insertAllSteps(self):
        self._initialize()
        for tomo in self.tomoList:
            tomoId = tomo.getTsId()
            self._insertFunctionStep(self._reconstructStep, tomoId)
            self._insertFunctionStep(self._createOutputStep, tomoId)

    # -------------------------- STEPS functions ------------------------------
    def _initialize(self):
        self.inReParticles = getattr(self.protPrepare.get(), prepareProtOutputs.relionParticles.name, None)
        self.outputSamplingRate = self.inReParticles.getTsSamplingRate() * self.binFactor.get()
        self.tomoSet = getObjFromRelation(self.inReParticles, self, SetOfTomograms)
        self.outTomoSet = SetOfTomograms.create(self._getPath(), template='tomograms%s.sqlite')
        self.outTomoSet.copyInfo(self.tomoSet)
        self.outTomoSet.setSamplingRate(self.outputSamplingRate)  # Bin factor provided is referred to the TS
        if self.recTomoMode.get() == SINGLE_TOMO:
            self.tomoList = [tomo.clone() for tomo in self.tomoSet if tomo.getTsId() == self.tomoId.get()]
        else:
            self.tomoList = [tomo.clone() for tomo in self.tomoSet]

    def _reconstructStep(self, tomoId):
        Plugin.runRelionTomo(self, 'relion_tomo_reconstruct_tomogram', self._genTomoRecCommand(tomoId))

    def _createOutputStep(self, tomoId):
        tomo = Tomogram()
        outFileName = self._getOutTomoFileName(tomoId)
        fixVolume(outFileName)
        tomo.setLocation(outFileName)
        tomo.setSamplingRate(self.outputSamplingRate)
        tomo.setOrigin()
        tomo.setTsId(tomoId)
        self.outTomoSet.append(tomo)
        self._defineOutputs(**{outputObjects.tomograms.name: self.outTomoSet})
        self._defineSourceRelation(self.inReParticles, self.outTomoSet)

    # -------------------------- INFO functions -------------------------------

    def _summary(self):
        summary = []
        if self.isFinished():
            if self.recTomoMode.get() == SINGLE_TOMO:
                summary.append('The selected tomogram was *%s*.' % self.tomoId.get())
        return summary

    # --------------------------- UTILS functions -----------------------------

    def _genTomoRecCommand(self, tomoId):
        cmd = '--t %s ' % self.inReParticles.getTomogramsStar()
        cmd += '--tn %s ' % tomoId
        cmd += '--o %s ' % self._getOutTomoFileName(tomoId)
        cmd += '--bin %.1f ' % self.binFactor.get()
        cmd += '--w %i ' % self.width.get()
        cmd += '--h %i ' % self.height.get()
        cmd += '--d %i ' % self.thickness.get()
        cmd += '--j %i ' % self.numberOfThreads.get()
        return cmd

    def _getOutTomoFileName(self, tomoId):
        return self._getExtraPath(tomoId + '.mrc')
















