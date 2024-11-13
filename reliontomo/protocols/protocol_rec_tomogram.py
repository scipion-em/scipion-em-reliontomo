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
from pyworkflow.protocol.params import FloatParam, IntParam, StringParam, PointerParam, EnumParam
from pyworkflow.utils import Message
from reliontomo import Plugin
from reliontomo.objects import RelionPSubtomogram
from reliontomo.protocols.protocol_base_relion import IS_RELION_50
from reliontomo.protocols.protocol_prepare_data import outputObjects as prepareProtOutputs
from reliontomo.utils import getProgram
from tomo.objects import Tomogram, SetOfTomograms
from tomo.utils import getObjFromRelation

# Reconstruct options
SINGLE_TOMO = 0
ALL_TOMOS = 1


class OutputObjects(Enum):
    tomograms = SetOfTomograms


class ProtRelionTomoReconstruct(EMProtocol):
    """ This protocol reconstructs a single tomogram using Relion. It is very useful
    to check if the protocol "Prepare data" has been applied correctly (in terms of flip
    options, for example).
    """
    _label = 'Reconstruct tomograms from prepare data prot'
    _possibleOutputs = OutputObjects

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
                      help='It is very useful to check if the protocol "Prepare data" has been applied correctly (in '
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

        form.addParam('binThreads', IntParam,
                      label='Relion threads',
                      default=3,
                      help='Number of threads used by Relion each time it is called in the protocol execution. For '
                           'example, if 2 Scipion threads and 3 Relion threads are set, the tomograms will be '
                           'processed in groups of 2 at the same time with a call of tomo3d with 3 threads each, so '
                           '6 threads will be used at the same time. Beware the memory of your machine has '
                           'memory enough to load together the number of tomograms specified by Scipion threads.')
        form.addParallelSection(threads=1, mpi=1)

    # -------------------------- INSERT steps functions -----------------------
    def _insertAllSteps(self):
        self._initialize()
        for tomo in self.tomoList:
            tomoId = tomo.getTsId()
            self._insertFunctionStep(self._reconstructStep, tomoId, needsGPU=False)
            self._insertFunctionStep(self._createOutputStep, tomoId, needsGPU=False)

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
            self.tomoList = [tomo.clone() for tomo in self.tomoSet if tomo.getTsId() in
                             self.inReParticles.getUniqueValues(RelionPSubtomogram.TS_ID_ATTRIBUTE)]

    def _reconstructStep(self, tomoId):
        nMpi = self.numberOfMpi.get()
        Plugin.runRelionTomo(self,
                             getProgram('relion_tomo_reconstruct_tomogram', nMpi=nMpi),
                             self._genTomoRecCommand(tomoId),
                             numberOfMpi=nMpi)

    def _createOutputStep(self, tomoId):
        tomo = Tomogram()
        outFileName = self._getOutTomoFileName(tomoId)
        fixVolume(outFileName)
        tomo.setLocation(outFileName)
        tomo.setSamplingRate(self.outputSamplingRate)
        tomo.setOrigin()
        tomo.setTsId(tomoId)
        self.outTomoSet.append(tomo)
        self._defineOutputs(**{OutputObjects.tomograms.name: self.outTomoSet})
        self._defineSourceRelation(self.inReParticles, self.outTomoSet)

    # -------------------------- INFO functions -------------------------------

    def _summary(self):
        summary = []
        if self.isFinished():
            if self.recTomoMode.get() == SINGLE_TOMO:
                summary.append('The selected tomogram was *%s*.' % self.tomoId.get())
        return summary

    # --------------------------- UTILS functions -----------------------------
    @classmethod
    def isDisabled(cls):
        """ Return True if this Protocol is disabled.
        Disabled protocols will not be offered in the available protocols."""
        return True if IS_RELION_50 else False

    def _genTomoRecCommand(self, tomoId):
        cmd = '--t %s ' % self.inReParticles.getTomogramsStar()
        cmd += '--tn %s ' % tomoId
        cmd += '--o %s ' % self._getOutTomoFileName(tomoId)
        cmd += '--bin %.1f ' % self.binFactor.get()
        cmd += '--w %i ' % self.width.get()
        cmd += '--h %i ' % self.height.get()
        cmd += '--d %i ' % self.thickness.get()
        cmd += '--j %i ' % self.binThreads.get()
        return cmd

    def _getOutTomoFileName(self, tomoId):
        return self._getExtraPath(tomoId + '.mrc')
















