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
import os
from enum import Enum
from os import mkdir
from os.path import join, exists
from imod.utils import generateDefocusIMODFileFromObject
from pwem.protocols import EMProtocol
from pyworkflow import BETA
from pyworkflow.object import Float
from pyworkflow.protocol import PointerParam, BooleanParam, LEVEL_ADVANCED
from reliontomo import Plugin
from reliontomo.constants import (IN_TOMOS_STAR, OUT_TOMOS_STAR, IN_COORDS_STAR,
                                  OPTIMISATION_SET_STAR)
from reliontomo.convert import writeSetOfTomograms, writeSetOfCoordinates
from reliontomo.objects import relionTomoMetadata

# data source
from tomo.utils import recoverTSFromObj

# Other constants
DEFOCUS = 'defocus'


class outputObjects(Enum):
    outputRelionParticles = relionTomoMetadata


class ProtRelionPrepareData(EMProtocol):
    """Prepare data for Relion 4
    """
    _label = 'Prepare data for Relion 4'
    _devStatus = BETA

    def __init__(self, **args):
        EMProtocol.__init__(self, **args)
        self.tsSet = None
        self.tomoSet = None
        self.coordScale = Float(1)

    # -------------------------- DEFINE param functions -----------------------

    def _defineParams(self, form):
        form.addSection(label='Input')

        form.addParam('inputCtfTs', PointerParam,
                      pointerClass='SetOfCTFTomoSeries',
                      label="Input CTF tomo series",
                      important=True,
                      allowsNull=False)
        form.addParam('handeness', BooleanParam,
                      label='Does focus decrease with Z distance?',
                      default=True,
                      expertLevel=LEVEL_ADVANCED,
                      help='It is the handedness of the tilt geometry and it is used to describe '
                           'whether the focus increases or decreases as a function of Z distance.'
                      )
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

        form.addParam('flipYZ', BooleanParam,
                       label='Has tomogram been flipped along Y and Z?',
                       default=False,
                       help='If the tomogram has been flipped along Y and Z (i.e. rotated around X) '
                            'after the reconstruction and before the particles have been picked, this '
                            'will apply the same transformation to the relion coordinate system. This will '
                            'allow relion to use particle positions defined in the X-rotated tomogram unchanged.')
        form.addParam('flipZ', BooleanParam,
                       label='Has the Z axis been flipped?',
                       default=False,
                       help='Same as above, in case the Z axis has been flipped. This can be used together with '
                            'the flipYZ option.')

    # -------------------------- INSERT steps functions -----------------------
    def _insertAllSteps(self):
        self._initialize()
        self._insertFunctionStep(self.convertInputStep)
        self._insertFunctionStep(self.relionImportTomograms)
        self._insertFunctionStep(self.relionImportParticles)
        self._insertFunctionStep(self.createOutputStep)

    # -------------------------- STEPS functions ------------------------------
    def _initialize(self):
        defocusDir = self._getExtraPath(DEFOCUS)
        if not exists(defocusDir):  # It can exist in case of Continue execution
            mkdir(defocusDir)
        self.coords = self.inputCoords.get()
        self.tsSet = recoverTSFromObj(self.coords, self)
        self.tomoSet = self.coords.getPrecedents()
        self.inputCtfs = self.inputCtfTs.get()
        # If coordinates are referred to a set of tomograms, they'll be rescaled
        # to be expressed in bin 1, as the ts images
        if self.tomoSet:
            self.coordScale.set(self.tomoSet.getSamplingRate() / self.tsSet.getSamplingRate())

    def convertInputStep(self):
        # Generate defocus files
        for ctfTomo in self.inputCtfs:
            defocusPath = self._getExtraPath(DEFOCUS, ctfTomo.getTsId())
            if not exists(defocusPath):  # It can exist in case of mode Continue execution
                mkdir(defocusPath)
            generateDefocusIMODFileFromObject(ctfTomo,
                                              join(defocusPath, ctfTomo.getTsId() + '.' + DEFOCUS),
                                              isRelion=True)
        thickness = self.coords.getPrecedents().getDim()[2]
        # Simulate the etomo files that serve as entry point to relion4
        self._simulateETomoFiles(self.tsSet, thickness=thickness, binned=self.coordScale,
                                 binByFactor=self.coordScale)
        # Write the tomograms star file
        writeSetOfTomograms(self.tsSet,
                            self._getStarFilename(IN_TOMOS_STAR),
                            prot=self,
                            ctfPlotterParentDir=self._getExtraPath(DEFOCUS),
                            eTomoParentDir=self._getTmpPath())
        # Write the particles star file
        writeSetOfCoordinates(self.inputCoords.get(),
                              self._getStarFilename(IN_COORDS_STAR),
                              sRate=self.tsSet.getSamplingRate(),
                              coordsScale=self.coordScale.get())

    def relionImportTomograms(self):
        Plugin.runRelionTomo(self, 'relion_tomo_import_tomograms', self._genImportTomosCmd())

    def relionImportParticles(self):
        Plugin.runRelionTomo(self, 'relion_tomo_import_particles', self._genImportSubtomosCmd())

    def createOutputStep(self):
        relionParticles = relionTomoMetadata(optimSetStar=self._getExtraPath(OPTIMISATION_SET_STAR),
                                             tsSamplingRate=self.tsSet.getSamplingRate(),
                                             relionBinning=self.coordScale.get(),
                                             nParticles=self.coords.getSize())

        self._defineOutputs(**{outputObjects.outputRelionParticles.name: relionParticles})
        self._store()

    # -------------------------- INFO functions -------------------------------
    def _validate(self):
        # TODO: generar los nombres culled --> tsId_culled.st:mrc cuando se quiten vistas con IMOD
        errorMsg = []
        tsSet = recoverTSFromObj(self.inputCoords.get(), self)
        if tsSet is None:
            errorMsg.append("Could not find any SetOfTiltSeries associated "
                            "with a transformation matrix")
        return errorMsg

    def _summary(self):
        msg = []
        if self.isFinished():
            if self.coordScale.get():
                msg.append('Coordinates were scaled using an scale factor of *%.2f* to be expressed considering the '
                           'size of the introduced tilt series' % self.coordScale.get())
        return msg

    # --------------------------- UTILS functions -----------------------------
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
        cmd = '--i %s ' % self._getStarFilename(IN_COORDS_STAR)
        cmd += '--o %s ' % self._getExtraPath()
        cmd += '--t %s ' % self._getStarFilename(OUT_TOMOS_STAR)
        if self.flipZCoords.get():
            cmd += '--flipZ '
        return cmd

    def _getStarFilename(self, fName):
        return self._getExtraPath(fName)

    def _decodeHandeness(self):
        return -1 if self.handeness.get() else 1

    def _simulateETomoFiles(self, imgSet, **kwargs):
        """Simulate the etomo files that serve as entry point to relion4
        """
        for ts in imgSet:
            # creating a folder where all data will be generate
            folderName = self._getTmpPath(ts.getTsId())
            makePath(folderName)
            # Create a symbolic link to the tomogram
            os.symlink(os.path.abspath(ts.getFirstItem().getFileName()),
                       os.path.join(folderName, ts.getTsId() + '.st'))
            ts.writeImodFiles(folderName, **kwargs)
