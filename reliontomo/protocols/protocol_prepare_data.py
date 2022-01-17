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
from os import mkdir, listdir
from os.path import join, exists
from imod.utils import generateDefocusIMODFileFromObject
from pwem.protocols import EMProtocol
from pyworkflow import BETA
from pyworkflow.protocol import PointerParam, PathParam, BooleanParam, LEVEL_ADVANCED, EnumParam
from reliontomo import Plugin
from reliontomo.constants import IN_TOMOS_STAR, IN_SUBTOMOS_STAR, OUT_TOMOS_STAR
from reliontomo.convert import writeSetOfTomograms, writeSetOfSubtomograms

# eTomo data source choices
ETOMO_FROM_PROT = 0
ETOMO_FROM_DIR = 1

# Other constants
DEFOCUS = 'defocus'


class ProtRelionPrepareData(EMProtocol):
    """Prepare data for Relion 4
    """
    _label = 'Prepare data for Relion 4'
    _devStatus = BETA

    def __init__(self, **args):
        EMProtocol.__init__(self, **args)
        self.tsSet = None
        self.tomoSet = None

    # -------------------------- DEFINE param functions -----------------------

    def _defineParams(self, form):
        form.addSection(label='CTFTomoSeries')

        form.addParam('inputCtfTs', PointerParam,
                      pointerClass='SetOfCTFTomoSeries',
                      label="Input CTF tomo series",
                      important=True,
                      allowsNull=False)
        form.addParam('handeness', BooleanParam,
                      label='Does focues decrease with Z distance?',
                      default=True,
                      expertLevel=LEVEL_ADVANCED,
                      help='It is the handedness of the tilt geometry and it is used to describe '
                           'whether the focus increases or decreases as a function of Z distance.'
                      )
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

        form.addSection(label='Coordinates')
        form.addParam('inputCoords', PointerParam,
                      pointerClass='SetOfCoordinates3D',
                      label="Input coordinates",
                      important=True,
                      allowsNull=False)
        form.addParam('flipZCoords', BooleanParam,
                      label='Flip Z coordinate?',
                      default=False,
                      expertLevel=LEVEL_ADVANCED
                      )

    # -------------------------- INSERT steps functions -----------------------
    def _insertAllSteps(self):
        self._initialize()
        self._insertFunctionStep(self.convertInputStep)
        self._insertFunctionStep(self.relionImportTomograms)
        self._insertFunctionStep(self.relionImportParticles)

    # -------------------------- STEPS functions ------------------------------
    def _initialize(self):
        defocusDir = self._getExtraPath(DEFOCUS)
        if not exists(defocusDir):  # It can exist in case of Continue execution
            mkdir(defocusDir)
        self.tsSet = self.inputCtfTs.get().getSetOfTiltSeries()
        self.tomoSet = self.inputCoords.get().getPrecedents()

    def convertInputStep(self):
        # Generate defocus files
        for ctfTomo in self.inputCtfTs.get():
            defocusPath = self._getExtraPath(DEFOCUS, ctfTomo.getTsId())
            if not exists(defocusPath):  # It can exist in case of mode Continue execution
                mkdir(defocusPath)
            generateDefocusIMODFileFromObject(ctfTomo,
                                              join(defocusPath, ctfTomo.getTsId() + '.' + DEFOCUS),
                                              isRelion=True)
        # Write the tomograms star file
        writeSetOfTomograms(self.tomoSet,
                            self._getStarFilename(IN_TOMOS_STAR),
                            prot=self,
                            tsSet=self.tsSet,
                            ctfPlotterParentDir=self._getExtraPath(DEFOCUS),
                            eTomoParentDir=self._getEtomoParentDir())
        # Write the particles star file
        writeSetOfSubtomograms(self.inputCoords.get(),
                               self._getStarFilename(IN_SUBTOMOS_STAR),
                               coordsScale=self.tomoSet.getSamplingRate() / self.tsSet.getSamplingRate())

    def relionImportTomograms(self):
        Plugin.runRelionTomo(self, 'relion_tomo_import_tomograms', self._genImportTomosCmd())

    def relionImportParticles(self):
        Plugin.runRelionTomo(self, 'relion_tomo_import_particles', self._genImportSubtomosCmd())

    # -------------------------- INFO functions -------------------------------
    def _validate(self):
        # TODO: generar los nombres culled --> tsId_culled.st:mrc cuando se quiten vistas con IMOD
        errorMsg = []
        # Check if files tilt.com and newst.com are contained in the eTomo corresponding directories
        NEWST_COM = 'newst.com'
        TILT_COM = 'tilt.com'
        tsIdList = [ts.getTsId() for ts in self.inputCtfTs.get().getSetOfTiltSeries()]
        eTomoParentDir = self._getEtomoParentDir()
        for tsId in tsIdList:
            currentTsEtomoDir = join(eTomoParentDir, tsId)
            currentTsErrMsg = ['\n%s [%s]:\n' % (tsId, currentTsEtomoDir)]
            if exists(currentTsEtomoDir):
                filesList = listdir(currentTsEtomoDir)
                if NEWST_COM not in filesList:
                    currentTsErrMsg.append('\t- File %s not found\n' % NEWST_COM)
                if TILT_COM not in filesList:
                    currentTsErrMsg.append('\t- File %s not found\n' % TILT_COM)
            else:
                currentTsErrMsg.append('\t- Directory %s not found\n' % currentTsEtomoDir)

            if len(currentTsErrMsg) > 1:
                errorMsg.append(''.join(currentTsErrMsg))

        return errorMsg

    def _summary(self):
        summary = []
        if self.isFinished():
            sRate = self.inputCoords.get().getPrecedents().getSamplingRate() / \
                    self.inputCtfTs.get().getSetOfTiltSeries().getSamplingRate()
            summary.append('Coordinates scaled to factor TS_SamplingRate / Coords_SamplingRate = %.2f' % sRate)

        return summary

    # --------------------------- UTILS functions -----------------------------
    def _getEtomoParentDir(self):
        if self.eTomoDataFrom.get() == ETOMO_FROM_PROT:
            return self.eTomoProt.get()._getExtraPath()
        else:
            return self.eTomoFilesPath.get()

    def _genImportTomosCmd(self):
        acq = self.tsSet.getAcquisition()
        cmd = '--i %s ' % self._getStarFilename(IN_TOMOS_STAR)
        cmd += '--o %s ' % self._getStarFilename(OUT_TOMOS_STAR)
        cmd += '--hand %s ' % self._decodeHandeness()
        cmd += '--angpix %s ' % self.tsSet.getSamplingRate()
        cmd += '--voltage %s ' % acq.getVoltage()
        cmd += '--Cs %s ' % acq.getSphericalAberration()
        cmd += '--Q0 %s ' % acq.getAmplitudeContrast()
        if self.flipYZ.get():
            cmd += '--flipYZ '
        if self.flipZ.get():
            cmd += '--flipZ '

        return cmd

    def _genImportSubtomosCmd(self):
        cmd = '--i %s ' % self._getStarFilename(IN_SUBTOMOS_STAR)
        cmd += '--o %s ' % self._getExtraPath()
        cmd += '--t %s ' % self._getStarFilename(OUT_TOMOS_STAR)
        if self.flipZCoords.get():
            cmd += '--flipZ '
        return cmd

    def _getStarFilename(self, fName):
        return self._getExtraPath(fName)

    def _decodeHandeness(self):
        return -1 if self.handeness.get() else 1
