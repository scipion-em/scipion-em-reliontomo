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
from emtable import Table

from pwem import EMObject
from pwem.objects import Volume, SetOfVolumes
from pyworkflow.object import String, Integer, Float
from relion.convert import OpticsGroups
from reliontomo.constants import OPT_TOMOS_STAR, OPT_PARTICLES_STAR, OPT_TRAJECTORIES_STAR, OPT_MANIFOLDS_STAR, \
    OPT_FSC_STAR
from tomo.objects import SubTomogram, SetOfSubTomograms
from tomo.protocols import ProtTomoBase


class RelionParticles(EMObject):

    def __init__(self, optimSetStar=None, samplingRate=None, nParticles=0, **kwargs):
        super().__init__(**kwargs)
        self._filesMaster = None
        self._nParticles = Integer(nParticles)
        self._tomograms = String()
        self._particles = String()
        self._trajectories = String()
        self._manifolds = String()
        self._referenceFsc = String()
        self._samplingRate = Float(samplingRate)
        if optimSetStar:
            self.filesMaster = optimSetStar

    def __str__(self):
        strRep = ''
        if self.getNumParticles():
            strRep += '%i items, ' % self.getNumParticles()
        if self.getSamplingRate():
            strRep += '%.2f Å/px' % self.getSamplingRate()
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

    def getSamplingRate(self):
        return self._samplingRate.get()

    def getNumParticles(self):
        return self._nParticles.get()

    def _readOptimSetStar(self, optimSetStar):
        dataTable = Table()
        dataTable.read(optimSetStar)
        rowObj = dataTable[0]  # This file only contains the different filenames related to the current STA step
        self._tomograms.set(rowObj.get(OPT_TOMOS_STAR, None))
        self._particles.set(rowObj.get(OPT_PARTICLES_STAR, None))
        self._trajectories.set(rowObj.get(OPT_TRAJECTORIES_STAR, None))
        self._manifolds.set(rowObj.get(OPT_MANIFOLDS_STAR, None))
        self._referenceFsc.set(rowObj.get(OPT_FSC_STAR, None))


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
    pass


# class PseudoSubtomogram(SubTomogram):
#
#     def __init__(self, tomoParticleName=None, randomSubset=None, opticsGroupId=None, xShiftAngst=None,
#                  yShiftAngst=None, zShiftAngst=None, ctfImage=None, **kwargs):
#         SubTomogram.__init__(self, **kwargs)
#         self._tomoParticleName = String(tomoParticleName)
#         self._randomSubset = Integer(randomSubset)
#         self._opticsGroupId = Integer(opticsGroupId)
#         self._xShiftAngst = Float(xShiftAngst)
#         self._yShiftAngst = Float(yShiftAngst)
#         self._zShiftAngst = Float(zShiftAngst)
#         self._ctfImage = String(ctfImage)
#
#     def getTomoParticleName(self):
#         return self._tomoParticleName.get()
#
#     def getRandomSubset(self):
#         return self._randomSubset.get()
#
#     def getOpticsGroupId(self):
#         return self._opticsGroupId.get()
#
#     def getShiftXAngst(self):
#         return self._xShiftAngst.get()
#
#     def getShiftYAngst(self):
#         return self._yShiftAngst.get()
#
#     def getShiftZAngst(self):
#         return self._zShiftAngst.get()
#
#     def getCtfImage(self):
#         return self._ctfImage.get()
#
#     def setTomoParticleName(self, val):
#         self._tomoParticleName.set(val)
#
#     def setRandomSubset(self, val):
#         self._randomSubset.set(val)
#
#     def setOpticsGroupId(self, val):
#         self._opticsGroupId.set(val)
#
#     def setShiftXAngst(self, val):
#         self._xShiftAngst.set(val)
#
#     def setShiftYAngst(self, val):
#         self._yShiftAngst.set(val)
#
#     def setShiftZAngst(self, val):
#         self._zShiftAngst.set(val)
#
#     def setCtfImage(self, val):
#         self._ctfImage.set(val)
#
#
# class SetOfPseudoSubtomograms(SetOfSubTomograms, ProtTomoBase):
#
#     def __init__(self, particlesStarFile=None, **kwargs):
#         SetOfSubTomograms.__init__(self, **kwargs)
#         self._starFile = String(particlesStarFile)
#
#     def getStarFile(self):
#         return self._starFile.get()
#
#     def setStarFile(self, starFile):
#         self._starFile.set(starFile)
#
#     def setStarFileAndOptics(self, starFile):
#         self._starFile.set(starFile)
#         self.getAcquisition().opticsGroupInfo.set(OpticsGroups.fromStar(starFile).toString())
#
#     def getOpticsGroupStr(self):
#         return self.getAcquisition().opticsGroupInfo.get()
#
#     def setOpticsGroupStr(self, particlesStarFile):
#         self.getAcquisition().opticsGroupInfo.set(OpticsGroups.fromStar(particlesStarFile).toString())
#
#     def copyInfo(self, other):
#         """ Copy basic information (sampling rate and ctf)
#         from other set of images to current one"""
#         super().copyInfo(other)
#         if hasattr(other, '_starFile'):
#             self.copyAttributes(other, '_starFile')
