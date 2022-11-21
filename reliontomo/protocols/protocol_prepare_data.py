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
from pyworkflow.utils import makePath, Message
from reliontomo import Plugin
from reliontomo.objects import createSetOfRelionPSubtomograms, RelionSetOfPseudoSubtomograms
from reliontomo.constants import (IN_TOMOS_STAR, OUT_TOMOS_STAR, IN_COORDS_STAR,
                                  OPTIMISATION_SET_STAR, OUT_PARTICLES_STAR, PSUBTOMOS_SQLITE)
from reliontomo.convert import writeSetOfTomograms, writeSetOfCoordinates, readSetOfPseudoSubtomograms
from reliontomo.utils import generateProjections
from tomo.utils import getNonInterpolatedTsFromRelations
import tomo.objects as tomoObj
from tomo.protocols.protocol_base import ProtTomoBase

# Other constants
DEFOCUS = 'defocus'
THICKNESS = 'thickness'
X_SIZE = 'x'
Y_SIZE = 'y'
DIMS = 'dims'


class outputObjects(Enum):
    relionParticles = RelionSetOfPseudoSubtomograms
    projected2DCoordinates = tomoObj.SetOfLandmarkModels


class ProtRelionPrepareData(EMProtocol, ProtTomoBase):
    """Prepare data for Relion 4
    """
    _label = 'Prepare data for Relion 4'
    _devStatus = BETA
    _possibleOutputs = outputObjects

    def __init__(self, **args):
        EMProtocol.__init__(self, **args)
        self.tsSet = None
        self.tomoSet = None
        self.tsReducedSet = None  # Reduced list of tiltseries with only those in the coordinates
        self.coordsVolIds = None  # Unique tomogram identifiers volId used in the coordinate set. In case is a subset
        self.coordScale = Float(1)

    # -------------------------- DEFINE param functions -----------------------

    def _defineParams(self, form):
        form.addSection(label=Message.LABEL_INPUT)

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
                           'whether the focus increases or decreases as a function of Z distance.')
        form.addParam('inputCoords', PointerParam,
                      pointerClass='SetOfCoordinates3D',
                      label="Input coordinates",
                      important=True,
                      allowsNull=False)
        form.addParam('inputTS', PointerParam,
                      help="Tilt series with alignment (non interpolated) used in the tomograms recnstruction. "
                           "To be deprecated!!",
                      pointerClass='SetOfTiltSeries',
                      label="Input tilt series",
                      important=True,
                      expertLevel=LEVEL_ADVANCED,
                      allowsNull=True)

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

        form.addParam('swapXY', BooleanParam,
                      label='Swap X with Y dimensions of the tilt series',
                      default=True,
                      help='This may be a trial and error parameter. Depending of the reconstruction path of '
                           'your tomograms we you may need to deactivate this option to get good results.')

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
        self.coordsVolIds = self._getTSIDFromCoordinates()
        self.tsSet = self._getTiltSeriesNonInterpolated()
        self.tomoSet = self.coords.getPrecedents()
        self.inputCtfs = self.inputCtfTs.get()
        # If coordinates are referred to a set of tomograms, they'll be rescaled
        # to be expressed in bin 1, as the ts images
        if self.tomoSet:
            self.coordScale.set(self.tomoSet.getSamplingRate() / self.tsSet.getSamplingRate())

    def _getTiltSeriesNonInterpolated(self):
        return self.inputTS.get() if self.inputTS.get() is not None else \
            getNonInterpolatedTsFromRelations(self.inputCoords.get(), self)

    def _getTSIDFromCoordinates(self):
        # TODO: Add this functionality to the SetOFCoordinates. May be useful for other cases
        volIds = self.coords.aggregate(["MAX"], "_tomoId", ["_tomoId"])
        volIds = [d['_tomoId'] for d in volIds]
        self.info("Involved tomograms in the coordinates: %s" % volIds)
        return volIds

    def convertInputStep(self):
        # Generate defocus files
        for ctfTomo in self.inputCtfs:
            defocusPath = self._getExtraPath(DEFOCUS, ctfTomo.getTsId())
            if not exists(defocusPath):  # It can exist in case of mode Continue execution
                mkdir(defocusPath)
            generateDefocusIMODFileFromObject(ctfTomo,
                                              join(defocusPath, ctfTomo.getTsId() + '.' + DEFOCUS),
                                              isRelion=True)
        # Thickness of the tomogram
        tomoSizeDict = {}
        tomoList = [tomo.clone() for tomo in self.coords.getPrecedents()]
        coordScale = self.coordScale.get()
        for tomo in tomoList:
            x, y, thickness = tomo.getDim()
            tomoSizeDict[tomo.getTsId()] = {X_SIZE: x * coordScale,
                                            Y_SIZE: y * coordScale,
                                            THICKNESS: thickness * coordScale}

        # Simulate the etomo files that serve as entry point to relion4
        self._simulateETomoFiles(self.tsSet, tomoSizeDict, binned=1, binByFactor=self.coordScale,
                                 whiteList=self.coordsVolIds, swapDims=self.swapXY.get())

        # Write the tomograms star file
        writeSetOfTomograms(self.tsSet,
                            self._getStarFilename(IN_TOMOS_STAR),
                            prot=self,
                            ctfPlotterParentDir=self._getExtraPath(DEFOCUS),
                            eTomoParentDir=self._getTmpPath(),
                            whiteList=self.coordsVolIds)
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
        # Pseudosubtomos
        coordSize = self.inputCoords.get().getBoxSize()

        psubtomoSet = createSetOfRelionPSubtomograms(self._getPath(),
                                                     self._getExtraPath(OPTIMISATION_SET_STAR),
                                                     template=PSUBTOMOS_SQLITE,
                                                     tsSamplingRate=self.tsSet.getSamplingRate(),
                                                     relionBinning=1,  # Coords are re-sampled to fit the TS size
                                                     boxSize=coordSize)
        # Fill the set with the generated particles
        readSetOfPseudoSubtomograms(psubtomoSet)

        self._defineOutputs(**{outputObjects.relionParticles.name: psubtomoSet})
        self._defineSourceRelation(self.inputCoords.get(), psubtomoSet)
        self._defineSourceRelation(self.inputCtfTs.get(), psubtomoSet)


        # Generate the fiducial model
        projections = generateProjections(self._getStarFilename(OUT_PARTICLES_STAR),
                                          self._getStarFilename(OUT_TOMOS_STAR))

        fiducialModelGaps = self._createSetOfLandmarkModels(suffix='Gaps')
        fiducialModelGaps.copyInfo(self.tsSet)
        fiducialModelGaps.setSetOfTiltSeries(self.tsSet)

        pos = 0
        for ts in self.tsSet:
            tsId = ts.getTsId()
            landmarkModelGapsFilePath = os.path.join(self._getExtraPath(),
                                                     str(tsId) + "_gaps.sfid")

            landmarkModelGaps = tomoObj.LandmarkModel(tsId=tsId,
                                              tiltSeriesPointer=ts,
                                              fileName=landmarkModelGapsFilePath,
                                              modelName=None,
                                              size=coordSize,
                                              applyTSTransformation=False)
            landmarkModelGaps.setTiltSeries(ts)

            while pos < len(projections) and projections[pos][0] == tsId:
                tiltIm = projections[pos][1] + 1
                chainId = projections[pos][2] + 1
                xCoor = int(round(projections[pos][3]))
                yCoor = int(round(projections[pos][4]))
                landmarkModelGaps.addLandmark(xCoor, yCoor, tiltIm,
                                              chainId, 0, 0)
                pos += 1
            fiducialModelGaps.append(landmarkModelGaps)

        self._defineOutputs(**{outputObjects.projected2DCoordinates.name: fiducialModelGaps})
        self._defineSourceRelation(self.tsSet,  fiducialModelGaps)
        self._store()

    # -------------------------- INFO functions -------------------------------
    def _validate(self):
        # TODO: generar los nombres culled --> tsId_culled.st:mrc cuando se quiten vistas con IMOD
        valMsg = []

        if not self.inputTS.get():
            try:
                self._getTiltSeriesNonInterpolated()
            except:
                valMsg.append('Unable to go via relations from the introduced coordinates to the '
                              'corresponding non-interpolated tilt series. Please introduce them using the '
                              'advanced parameter "Tilt series with alignment..."')
        return valMsg

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

    def _simulateETomoFiles(self, imgSet, tomoSizeDict, whiteList=None, **kwargs):
        """Simulate the etomo files that serve as entry point to relion4
        """
        for ts in imgSet:
            tsId = ts.getTsId()
            # Get the size of the corresponding tomogram, as it may differ from one to another
            tomoIdMatchDict = tomoSizeDict.get(tsId, None)
            if tomoIdMatchDict and (whiteList is None or tsId in whiteList):
                kwargs[THICKNESS] = tomoIdMatchDict[THICKNESS]
                kwargs[DIMS] = (tomoIdMatchDict[X_SIZE], tomoIdMatchDict[Y_SIZE])
                # creating a folder where all data will be generate
                folderName = self._getTmpPath(tsId)
                makePath(folderName)
                # Create a symbolic link to the tiltseries image file
                os.symlink(os.path.abspath(ts.getFirstItem().getFileName()),
                           os.path.join(folderName, ts.getTsId() + '.st'))
                ts.writeImodFiles(folderName, **kwargs)
