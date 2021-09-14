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


class PseudoSubtomogram(SubTomogram):

    def __init__(self, **kwargs):
        SubTomogram.__init__(self, **kwargs)
        self._tomoParticleName = String()
        self._randomSubset = Integer()
        self._opticsGroupId = Integer()
        self._xShiftAngst = Float()
        self._yShiftAngst = Float()
        self._zShiftAngst = Float()
        self._ctfImage = String()

    # TODO: add all the getters and setters


class SetOfPseudoSubtomograms(SetOfSubTomograms):
    def __init__(self, particlesStarFile, **kwargs):
        SetOfSubTomograms.__init__(self, **kwargs)
        self._opticsGroup = OpticsGroups.fromStar(particlesStarFile)

    # TODO: add all the getters and setters
