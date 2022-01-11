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
from os import mkdir

from pyworkflow.utils import getParentFolder
from reliontomo.protocols.protocol_base_import_from_star import ProtBaseImportFromStar
from tomo.objects import SetOfSubTomograms


class ProtImportSubtomogramsFromStar(ProtBaseImportFromStar):
    """Protocol to import a set of subtomograms from a star file"""

    _label = 'import subtomograms from a star file'

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.linkedStarFileName = 'inSubtomograms.star'
        self.linkedSubtomosDirName = 'subtomograms'
        self.isSubtomoStarFile = True

    # --------------------------- STEPS functions -----------------------------

    def _initialize(self):
        super()._initialize()
        # Generate the firectoy in which the linked subtomograms pointed from the star file will be stored
        mkdir(self._getExtraPath(self.linkedSubtomosDirName))

    def _importStep(self):
        # Generate the corresponding precedents and 3d coordinates
        super()._importStep()
        # Generate the set of subtomograms
        outputTomoSet = getattr(self, 'outputTomograms', None)
        outputCoordsSet = getattr(self, 'outputCoordinates', None)
        subtomoSet = SetOfSubTomograms.create(self._getPath(), template='setOfSubTomograms%s.sqlite')
        subtomoSet.setSamplingRate(self.sRate)
        subtomoSet.setAcquisition(outputTomoSet.getAcquisition())
        self.reader.starFile2Subtomograms(subtomoSet, outputCoordsSet, self._getExtraPath(self.linkedSubtomosDirName),
                                          getParentFolder(self.starFile.get()))
        self._defineOutputs(outputSubtomograms=subtomoSet)

