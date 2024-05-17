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
from reliontomo.protocols.protocol_base_import_from_star import ProtBaseImportFromStar
from tomo.objects import SetOfTomograms, SetOfCoordinates3D


class outputObjects(Enum):
    tomograms = SetOfTomograms()
    coordinates = SetOfCoordinates3D()


class ProtImportCoordinates3DFromStar(ProtBaseImportFromStar):
    """Protocol to import a 3D coordinates from a relion star file as the one provided in the tutorial"""

    _label = 'import 3D coordinates from a star file'
    _possibleOutputs = outputObjects

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.isCoordsFile = True
