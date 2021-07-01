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
from pyworkflow.protocol import PointerParam, PathParam, BooleanParam, LEVEL_ADVANCED, EnumParam
from reliontomo import Plugin
from reliontomo.constants import IN_TOMOS_STAR, IN_SUBTOMOS_STAR, OUT_TOMOS_STAR
from reliontomo.convert import writeSetOfTomograms, writeSetOfSubtomograms

# eTomo data source choices
ETOMO_FROM_PROT = 0
ETOMO_FROM_DIR = 1


class ProtRelionPrepareData(EMProtocol):
    """
    """
    _label = 'Load and prepare data as expected by relion 4'
    _devStatus = BETA
    tomoSet = None
    acquisition = None

    # -------------------------- DEFINE param functions -----------------------

    def _defineParams(self, form):
        form.addSection(label='CTFTomoSeries')

        form.addParam('inputCTFTS', PointerParam,
                      pointerClass='SetOfCTFTomoSeries',
                      label="Input CTF tomo series",
                      important=True,
                      allowsNull=False,
                      help='Select the input set of CTF tomo series from the project.')
        form.addParam('handeness', BooleanParam,
                      label='Does focues decrease with Z distance?',
                      default=True,
                      expertLevel=LEVEL_ADVANCED,
                      help='It is the handedness of the tilt geometry and it is used to describe '
                           'whether the focus increases or decreases as a function of Z distance.'
                      )
        # form.addParam('ctfPlotterFilesPath', PathParam,
        #               label="IMOD's CTFPLotter results parent directory",
        #               important=True,
        #               allowsNull=False,
        #               help='Used to provide the *.defocus* files obtained with CTFPlotter. '
        #                    'It is expected to be the parent folder where all your resulting subdirectories '
        #                    'obtained with IMOD-CTFPlotter are contained.\n\n'
        #                    'There *must* be one subdirectory '
        #                    '*per tilt series* and each file *must* contain as many rows as the '
        #                    'number of tilt images which compose the corresponding tilt series.\n\n'
        #                    '*Example:*\n'
        #                    'CTF manual estimation protocol folder:\n'
        #                    '    |_extra\n'
        #                    '        |_ts1\n'
        #                    '            |_ts1.defocus\n'
        #                    '        |_ts2\n'
        #                    '            |_ts2.defocus\n\n'
        #                    'In this case, the CTFPlotter directory would be the path to\n'
        #                    '*CTF manual estimation protocol folder*.'
        #               )

        group = form.addGroup('IMOD related arguments')
        group.addParam('eTomoDataFrom', EnumParam,
                       label='Choose IMOD-eTomo data source',
                       display=EnumParam.DISPLAY_HLIST,
                       choices=['Protocol', 'Directory'],
                       important=True,
                       default=ETOMO_FROM_PROT)
        group.addParam('eTomoProt', PointerParam,
                       pointerClass='ProtImodEtomo',
                       label='IMOD-eTomo protocol',
                       important=True,
                       condition='eTomoDataFrom == %s' % ETOMO_FROM_PROT)
        group.addParam('eTomoFilesPath', PathParam,
                       label="IMOD's eTomo results",
                       important=True,
                       allowsNull=False,
                       condition='eTomoDataFrom == %s' % ETOMO_FROM_DIR,
                       help='It is expected to be the parent folder where all your resulting subdirectories '
                            'obtained with IMOD-eTomo are contained.\n\n'
                            'There *must* be one subdirectory per *tilt series*, and each of one '
                            '*must* contain the corresponding *newst.com* and *tilt.com files*.\n\n'
                            '*Note:*\n'
                            '*Example:*\n'
                            'IMOD - eTomo manual estimation protocol folder:\n'
                            '    |_extra\n'
                            '        |_ts1\n'
                            '            |_newst.com\n'
                            '            |_tilt.com\n'  
                            '        |_ts2\n'
                            '            |_newst.com\n'
                            '            |_tilt.com\n' 
                            'In Scipion case, the IMOD - eTomo directory would be the path to\n'
                            '*[IMOD - eTomo protocol]/extra*.')
        group.addParam('flipYZ', BooleanParam,
                       label='Has IMOD tomogram been flipped along Y and Z?',
                       default=False,
                       help='If the IMOD tomogram has been flipped along Y and Z (i.e. rotated around X) '
                            'after an IMOD reconstruction and before the particles have been picked, this '
                            'will apply the same transformation to the relion coordinate system. This will '
                            'allow relion to use particle positions defined in the X-rotated tomogram unchanged.')
        group.addParam('flipZ', BooleanParam,
                       label='Has the Z axis been flipped?',
                       default=False,
                       help='Same as above, in case the Z axis has been flipped. This can be used together with '
                            'the flipYZ option.')

        form.addSection(label='Subtomograms')
        form.addParam('inputSubtomos', PointerParam,
                      pointerClass='SetOfSubTomograms',
                      label="Input set of subtomograms",
                      important=True,
                      allowsNull=False,
                      help='Select the input set of subtomograms from the project.')
        form.addParam('flipZCoords', BooleanParam,
                      label='Flip Z coordinate?',
                      default=False,
                      expertLevel=LEVEL_ADVANCED
                      )

    def _insertAllSteps(self):
        self._initialize()
        self._insertFunctionStep(self._convertInputStep)
        self._insertFunctionStep(self._relionImportTomograms)
        self._insertFunctionStep(self._relionImportParticles)

    def _initialize(self):
        self.tomoSet = self.inputSubtomos.get().getCoordinates3D().getPrecedents()
        self.acquisition = self.tomoSet.getAcquisition()

    def _convertInputStep(self):
        # Write the tomograms star file
        writeSetOfTomograms(self.tomoSet, self._getInTomosStarFilename(), prot=self, tsSet=self.inputTS.get(),
                            ctfPlotterDir=self.ctfPlotterFilesPath.get(), eTomoDir=self.eTomoFilesPath.get())
        # Write the particles star file
        writeSetOfSubtomograms(self.inputSubtomos.get(), self._getInSubtomosStarFilename())

    def _relionImportTomograms(self):
        self.runJob('relion_tomo_import_tomograms', self._genImportTomosCmd(), env=Plugin.getEnviron())

    def _relionImportParticles(self):
        self.runJob('flipZCoords', self._genImportSubtomosCmd(), env=Plugin.getEnviron())

    def _validate(self):
        # TODO: generar cada .defocus file a partir de la CTFTomoSeries --> Fede
        # TODO: lo mismo con los tilt.com y newst.com
        # TODO: el orderlist tiene que ir en un csv con order, ángulo
        # TODO: generar los nombres culled --> tsId_culled.st:mrc
        pass

    def _genImportTomosCmd(self):
        cmd = '--i %s ' % self._getStarFilename(IN_TOMOS_STAR)
        cmd += '--o %s ' % self._getStarFilename(OUT_TOMOS_STAR)
        cmd += '--hand %s ' % self._decodeHandeness()
        cmd += '--angpix %s ' % self.inputSubtomos.get().getSamplingRate()
        cmd += '--voltage %s ' % self.acquisition.getVoltage()
        cmd += '--Cs %s ' % self.acquisition.getSphericalAberration()
        cmd += '--Q0 %s ' % self.acquisition.getAmplitudeContrast()
        if self.flipYZ.get():
            cmd += '--flipYZ '
        if self.flipZ.get():
            cmd += '--flipZ '

        return cmd

    def _genImportSubtomosCmd(self):
        cmd = '--i %s ' % self._getStarFilename(IN_SUBTOMOS_STAR)
        cmd += '--t %s ' % self._getStarFilename(OUT_TOMOS_STAR)
        cmd += '--o %s ' % self._getExtraPath()
        if self.flipZCoords().get():
            cmd += '--flipZ '
        return cmd

    def _getStarFilename(self, fName):
        return self._getExtraPath(fName)

    def _decodeHandeness(self):
        return -1 if self.handeness.get() else 1
