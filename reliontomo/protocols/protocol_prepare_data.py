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
"""
This module contains the protocol for 3d classification with relion.
"""
from pwem.protocols import EMProtocol
from pyworkflow import BETA
from pyworkflow.protocol import PointerParam, PathParam
from reliontomo.constants import IN_TOMO_STAR, IN_SUBTOMOS_STAR
from reliontomo.convert import writeSetOfTomograms


class ProtRelionPrepareData(EMProtocol):
    """
    """
    _label = 'Load and prepare data as expected by relion 4'
    _devStatus = BETA
    # -------------------------- DEFINE param functions -----------------------

    def _defineParams(self, form):
        form.addSection(label='Input')

        form.addParam('inputTS', PointerParam,
                      pointerClass='SetOfTiltSeries',
                      label="Input tilt series",
                      important=True,
                      allowsNull=False,
                      help='Select the input tilt series from the project.')
        form.addParam('inputSubtomos', PointerParam,
                      pointerClass='SetOfSubTomograms',
                      label="Input set of subtomograms",
                      important=True,
                      allowsNull=False,
                      help='Select the input set of subtomograms from the project.')
        form.addParam('ctfPlotterFilesPath', PathParam,
                      label="IMOD's CTFPLotter results parent directory",
                      important=True,
                      allowsNull=False,
                      help='Used to provide the *.defocus* files obtained with CTFPlotter. '
                           'It is expected to be the parent folder where all your resulting subdirectories '
                           'obtained with IMOD-CTFPlotter are contained.\n\n'
                           'There *must* be one subdirectory '
                           '*per tilt series* and each file *must* contain as many rows as the '
                           'number of tilt images which compose the corresponding tilt series.\n\n'
                           '*Example:*\n'
                           'CTF manual estimation protocol folder:\n'
                           '    |_extra\n'
                           '        |_ts1\n'
                           '            |_ts1.defocus\n'
                           '        |_ts2\n'
                           '            |_ts2.defocus\n\n'
                           'In this case, the CTFPlotter directory would be the path to\n'
                           '*CTF manual estimation protocol folder*.'
                      )
        form.addParam('eTomoFilesPath', PathParam,
                      label="IMOD's eTomo results parent directory",
                      important=True,
                      allowsNull=False,
                      help='It is expected to be the parent folder where all your resulting subdirectories '
                           'obtained with IMOD-eTomo are contained.\n\n'
                           'There *must* be one subdirectory per *tilt series*, and each of one '
                           '*must* contain the corresponding *newst.com* and *tilt.com files*.\n\n'
                           '*Note:*\n'
                           'For a better understanding, check the help of the parameter CTFPlotter '
                           'results parent directory.')

    def _insertAllSteps(self):
        self._insertFunctionStep(self._convertInputStep)
        self._insertFunctionStep(self._relionImportTomograms)
        self._insertFunctionStep(self._relionImportParticles)
        self._insertFunctionStep(self._createOutputStep)

    def _convertInputStep(self):
        # Write the tomograms star file
        tomoSet = self.inputSubtomos.get().getCoordinates3D().getPrecedents()
        writeSetOfTomograms(tomoSet, self._getInTomosStarFilename(), prot=self, tsSet=self.inputTS.get(),
                            ctfPlotterDir=self.ctfPlotterFilesPath.get(), eTomoDir=self.eTomoFilesPath.get())
        # Write the particles star file
        # Writer.writeSetOfSubtomograms()

    def _relionImportTomograms(self):
        pass

    def _relionImportParticles(self):
        pass

    def _createOutputStep(self):
        pass

    def _validate(self):
        # TODO: asociar cada TS con cada .defocus file
        # TODO: lo mismo con los tilt.com y newst.com
        # TODO: el orderlist tiene que ir en un csv con order, ángulo
        # TODO: generar los nombres culled --> tsId_culled.st:mrc
        pass

    def _getInTomosStarFilename(self):
        return self._getExtraPath(IN_TOMO_STAR)

    def _getInSubtomosStarFilename(self):
        return self._getExtraPath(IN_SUBTOMOS_STAR)

