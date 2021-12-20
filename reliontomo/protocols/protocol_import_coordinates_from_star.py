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
from reliontomo.protocols.protocol_base_import_from_star import ProtBaseImportFromStar
from tomo.objects import SetOfTomograms, TomoAcquisition


class ProtImportCoordinates3DFromStar(ProtBaseImportFromStar):
    """Protocol to import a set of 3D coordinates from a star file"""

    _label = 'import coordinates 3D from a star file'

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.linkedStarFileName = 'in3dCoordinates.star'

    # --------------------------- STEPS functions -----------------------------

    # def _importStep(self):
    #     # Generate the precedents (set of tomograms which the coordinates are referred to) if necessary
    #     tomoPrecedentsSet = SetOfTomograms.create(self._getPath(), template='tomograms%s.sqlite')
    #     tomoPrecedentsSet.setSamplingRate(self.sRate)
    #     tomoPrecedentsSet.setAcquisition(TomoAcquisition(angleMin=-60, angleMax=60, step=3))  # Generic values
    #     self._fillPrecedentsSet(tomoPrecedentsSet)
    #     self._defineOutputs(outputSetOfTomograms=tomoPrecedentsSet)
    #
    #     # Read the star file and generate the corresponding set
    #     coordSet = self._createSetOfCoordinates3D(tomoPrecedentsSet)
    #     coordSet.setSamplingRate(self.sRate)
    #     coordSet.setBoxSize(self.boxSize.get())
    #     self.reader.starFile2Coords3D(coordSet, tomoPrecedentsSet)
    #
    #     self._defineOutputs(outputCoordinates=coordSet)
    #     self._defineSourceRelation(tomoPrecedentsSet, coordSet)