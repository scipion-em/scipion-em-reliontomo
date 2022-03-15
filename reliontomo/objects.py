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
from os.path import exists, join

from emtable import Table
from pwem import EMObject
from pwem.objects import Volume, SetOfVolumes
from pyworkflow.object import String, Integer, Float
from reliontomo.constants import OPT_TOMOS_STAR, OPT_PARTICLES_STAR, OPT_TRAJECTORIES_STAR, OPT_MANIFOLDS_STAR, \
    OPT_FSC_STAR, OUT_TOMOS_STAR, OUT_PARTICLES_STAR, TRAJECTORIES_STAR, MANIFOLDS_STAR, \
    FSC_REF_STAR


class EnumRe4GenFilesProps(Enum):
    _tomograms = OUT_TOMOS_STAR
    _particles = OUT_PARTICLES_STAR
    _trajectories = TRAJECTORIES_STAR
    _manifolds = MANIFOLDS_STAR
    _referenceFsc = FSC_REF_STAR


class relionTomoMetadata(EMObject):

    def __init__(self, optimSetStar=None, relionBinning=None, tsSamplingRate=None, nParticles=0, **kwargs):
        super().__init__(**kwargs)
        self._filesMaster = None
        self._nParticles = Integer(nParticles)
        self._tomograms = String()
        self._particles = String()
        self._trajectories = String()
        self._manifolds = String()
        self._referenceFsc = String()
        self._relionBinning = Float(relionBinning)
        self._tsSamplingRate = Float(tsSamplingRate)
        if optimSetStar:
            self.filesMaster = optimSetStar

    def __str__(self):
        strRep = ''
        if self.getNumParticles():
            strRep += '%i items, ' % self.getNumParticles()
        if self.getRelionBinning():
            strRep += 'binning %.1f' % self.getRelionBinning()
        # if self.getTsSamplingRate():
        #     strRep += 'TS Sampling Rate %.2f Å/px' % self.getTsSamplingRate()
        return self.getClassName() + ' %s' % ('(%s)' % strRep if strRep else '')

    @property
    def filesMaster(self):
        return self._filesMaster

    @filesMaster.setter
    def filesMaster(self, optimSetStar):
        try:
            self._readOptimSetStar(optimSetStar)
            self._filesMaster = optimSetStar
        except FileNotFoundError:
            raise FileNotFoundError('Unable to find file %s' % optimSetStar)

    def updateGenFiles(self, extraPath):
        """Some protocols don't generate the optimisation_set.star file. In that case, the input Object which
        represents it will be copied and, after that, this method will be used to update the corresponding
        attributes with the generated files."""
        for p in EnumRe4GenFilesProps:
            currentFile = join(extraPath, p.value)
            if exists(currentFile):
                setattr(self, p.name, String(currentFile))

    def getTomograms(self):
        return self._tomograms.get()

    def getParticles(self):
        return self._particles.get()

    def getTrajectories(self):
        return self._trajectories.get()

    def getManifolds(self):
        return self._manifolds.get()

    def getReferenceFsc(self):
        return self._referenceFsc.get()

    def getRelionBinning(self):
        return self._relionBinning.get()

    def getTsSamplingRate(self):
        return self._tsSamplingRate.get()

    def getNumParticles(self):
        return self._nParticles.get()

    def getCurrentSamplingRate(self):
        return self.getTsSamplingRate() * self.getRelionBinning()

    def _readOptimSetStar(self, optimSetStar):
        dataTable = Table()
        dataTable.read(optimSetStar)
        rowObj = dataTable[0]  # This file only contains the different filenames related to the current STA step
        self._tomograms.set(rowObj.get(OPT_TOMOS_STAR, None))
        self._particles.set(rowObj.get(OPT_PARTICLES_STAR, None))
        self._trajectories.set(rowObj.get(OPT_TRAJECTORIES_STAR, None))
        self._manifolds.set(rowObj.get(OPT_MANIFOLDS_STAR, None))
        self._referenceFsc.set(rowObj.get(OPT_FSC_STAR, None))

    def copyInfo(self, other):
        self.copyAttributes(other, '_nParticles', '_tomograms', '_particles', '_trajectories',
                            '_manifolds', '_referenceFsc', '_relionBinning', '_tsSamplingRate')


class PSubtomogram(Volume):

    def __init__(self, fileName=None, samplingRate=None, ctfFile=None, tsId=None, classId=None, **kwargs):
        super().__init__(**kwargs)
        if fileName:
            self.setFileName(fileName)
        if samplingRate:
            self.setSamplingRate(samplingRate)
        if classId:
            self.setClassId(classId)
        self.tsId = String(tsId)
        self.ctfFile = String(ctfFile)


class SetOfPseudoSubtomograms(SetOfVolumes):
    ITEM_TYPE = PSubtomogram
