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

from reliontomo.constants import IN_PARTICLES_STAR
from reliontomo.protocols.protocol_re5_base_import_from_star import ProtBaseRe5ImportFromStar
from tomo.objects import SetOfTomograms, SetOfCoordinates3D


class outputObjects(Enum):
    tomograms = SetOfTomograms()
    coordinates = SetOfCoordinates3D()


class ProtImportRe5Coordinates3DFromStar(ProtBaseRe5ImportFromStar):
    """Protocol to import a 3D coordinates from a relion star file as the one provided in the tutorial"""

    _label = 'import re5 3D coordinates'
    _possibleOutputs = outputObjects

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.linkedStarFileName = IN_PARTICLES_STAR
