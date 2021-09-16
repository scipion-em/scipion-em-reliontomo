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
from pyworkflow.object import String, Integer, Float
from relion.convert import OpticsGroups
from tomo.objects import SubTomogram, SetOfSubTomograms
from tomo.protocols import ProtTomoBase


class PseudoSubtomogram(SubTomogram):

    def __init__(self, tomoParticleName=None, randomSubset=None, opticsGroupId=None, xShiftAngst=None,
                 yShiftAngst=None, zShiftAngst=None, ctfImage=None, **kwargs):
        SubTomogram.__init__(self, **kwargs)
        self._tomoParticleName = String(tomoParticleName)
        self._randomSubset = Integer(randomSubset)
        self._opticsGroupId = Integer(opticsGroupId)
        self._xShiftAngst = Float(xShiftAngst)
        self._yShiftAngst = Float(yShiftAngst)
        self._zShiftAngst = Float(zShiftAngst)
        self._ctfImage = String(ctfImage)

    def getTomoParticleName(self):
        return self._tomoParticleName.get()

    def getRandomSubset(self):
        return self._randomSubset.get()

    def getOpticsGroupId(self):
        return self._opticsGroupId.get()

    def getShiftXAngst(self):
        return self._xShiftAngst.get()

    def getShiftYAngst(self):
        return self._yShiftAngst.get()

    def getShiftZAngst(self):
        return self._zShiftAngst.get()

    def getCtfImage(self):
        return self._ctfImage.get()

    def setTomoParticleName(self, val):
        self._tomoParticleName.set(val)

    def setRandomSubset(self, val):
        self._randomSubset.set(val)

    def setOpticsGroupId(self, val):
        self._opticsGroupId.set(val)

    def setShiftXAngst(self, val):
        self._xShiftAngst.set(val)

    def setShiftYAngst(self, val):
        self._yShiftAngst.set(val)

    def setShiftZAngst(self, val):
        self._zShiftAngst.set(val)

    def setCtfImage(self, val):
        self._ctfImage.set(val)


class SetOfPseudoSubtomograms(SetOfSubTomograms, ProtTomoBase):

    def __init__(self, particlesStarFile=None, **kwargs):
        SetOfSubTomograms.__init__(self, **kwargs)
        self._opticsGroup = OpticsGroups.fromStar(particlesStarFile) if particlesStarFile else None

    def getOpticsGroupObj(self):
        return self._opticsGroup.get()

    def setOpticsGroupObjFromStar(self, particlesStarFile):
        self._opticsGroup = OpticsGroups.fromStar(particlesStarFile)

    def setOpticsGroupObj(self, opticsGroupObj):
        self._opticsGroup = opticsGroupObj
