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
from pyworkflow.protocol import PointerParam, BooleanParam, LEVEL_ADVANCED, IntParam
from pyworkflow.utils import makePath, Message
from reliontomo import Plugin
from reliontomo.objects import createSetOfRelionPSubtomograms, RelionSetOfPseudoSubtomograms
from reliontomo.constants import (IN_TOMOS_STAR, OUT_TOMOS_STAR, IN_COORDS_STAR,
                                  OPTIMISATION_SET_STAR, OUT_PARTICLES_STAR, PSUBTOMOS_SQLITE, IN_PARTICLES_STAR)
from reliontomo.convert import writeSetOfTomograms, writeSetOfCoordinates, readSetOfPseudoSubtomograms, convert50_tomo
from reliontomo.protocols.protocol_re5_base_extract_subtomos_and_rec_particle import \
    ProtRelion5ExtractSubtomoAndRecParticleBase
from reliontomo.utils import generateProjections
import tomo.objects as tomoObj


# # Other constants
# DEFOCUS = 'defocus'
# THICKNESS = 'thickness'
# X_SIZE = 'x'
# Y_SIZE = 'y'
# DIMS = 'dims'
# SHIFTS = 'shift'


# class outputObjects(Enum):
#     relionParticles = RelionSetOfPseudoSubtomograms
#     projected2DCoordinates = tomoObj.SetOfLandmarkModels


class ProtRelion5ExtractSubtomos(ProtRelion5ExtractSubtomoAndRecParticleBase):
    """extracts the relevant cropped areas of the tilt series images for each individual particle and saves them as
    CTF-premultiplied extracted 2D image stacks (or as 3D volumes).
    """
    _label = 'Extract subtomos'
    _devStatus = BETA

    # _possibleOutputs = outputObjects

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        # self.tsSet = None
        # self.tomoSet = None
        # self.tsReducedSet = None  # Reduced list of tiltseries with only those in the coordinates
        # self.matchingTSIds = None  # Unique tomogram identifiers volId used in the coordinate set. In case is a subset
        self.coordScale = Float(1)
        self.tsDict = dict()
        self.ctfDict = dict()
        self.tomoDict = dict()

    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        form.addSection(label=Message.LABEL_INPUT)
        form.addParam('inputCoords', PointerParam,
                      pointerClass='SetOfCoordinates3D',
                      label="Coordinates",
                      important=True,
                      allowsNull=False)
        form.addParam('inputCtfTs', PointerParam,
                      pointerClass='SetOfCTFTomoSeries',
                      label="CTF tomo series",
                      important=True,
                      allowsNull=False)
        form.addParam('inputTS', PointerParam,
                      help="Tilt series with alignment (non interpolated) used in the tomograms reconstruction.",
                      pointerClass='SetOfTiltSeries',
                      label="Tilt series",
                      important=True,
                      allowsNull=False)
        form.addParam('handeness', BooleanParam,
                      label='Does focus decrease with Z distance?',
                      default=True,
                      help='It is the handedness of the tilt geometry and it is used to describe '
                           'whether the focus increases or decreases as a function of Z distance.')
        form.addSection(label='Reconstruct')
        self._defineCommonRecParams(form)
        form.addParam('maxDose', IntParam,
                      label='Maximum dose (e/Å^2)',
                      default=50,
                      help='Tilt series frames with a dose higher than this maximum dose (in electrons per squared '
                           'Angstroms) will not be included in the 3D pseudo-subtomogram, or in the 2D stack. For the '
                           'latter, this will affect disc I/O operations and increase speed.')
        form.addParam('minNoFrames', IntParam,
                      label='Minimum no. frames',
                      default=1,
                      help='Each selected pseudo-subtomogram need to be visible in at least this number of tilt series '
                           'frames with doses below the maximum dose.')
        form.addParam('write2dStacks', BooleanParam,
                      label='Write output as 2D stacks?',
                      default=True,
                      help='If set to Yes, this program will write output subtomograms as 2D substacks. This is new '
                           'as of relion-4.1, and the preferred way of generating subtomograms. If set to No, '
                           'then relion-4.0 3D pseudo-subtomograms will be written out. Either can be used in '
                           'subsequent refinements and classifications..')
        form.addParam('outputInFloat16', BooleanParam,
                      label='Write output in float16?',
                      expertLevel=LEVEL_ADVANCED,
                      default=False,
                      help='If set to Yes, this program will write output images in float16 MRC format. This will '
                           'save a factor of two in disk space compared to the default of writing in float32. '
                           'Note that RELION and CCPEM will read float16 images, but other programs may '
                           'not (yet) do so.')

    # -------------------------- INSERT steps functions -----------------------
    def _insertAllSteps(self):
        self._initialize()
        self._insertFunctionStep(self.convertInputStep)
        self._insertFunctionStep(self.extractSubtomos)
    #     self._insertFunctionStep(self.createOutputStep)
    #
    # -------------------------- STEPS functions ------------------------------
    def _initialize(self):
        coords = self.inputCoords.get()
        tsSet = self.inputTS.get()
        ctfSet = self.inputCtfTs.get()

        # The ccordinates need to be re-scaled to be at the same size of the tilt-series
        self.coordScale.set(coords.getSamplingRate() / tsSet.getSamplingRate())

        # Compute matching TS id among coordinates, the tilt-series and the CTFs, they all could be a subset
        coordsTsIds = coords.getTSIds()
        self.info("TsIds present in coordinates are: %s" % coordsTsIds)
        tsIds = tsSet.getTSIds()
        self.info("TsIds present in Tilt series are: %s" % tsIds)
        ctfTsIds = ctfSet.getTSIds()
        presentTsIds = set(coordsTsIds) & set(tsIds) & set(ctfTsIds)

        # Validate the intersection
        if len(presentTsIds) > 0:
            self.info("Tilt series associated in coordinates, CTFs, and tilt series are %s" % presentTsIds)
        else:
            raise Exception("There isn't any common tilt-series ids among the coordinates, CTFs, and tilt-series "
                            "introduced.")

        # Manage the TS, CTF tomo Series and tomograms
        self.tsDict = {ts.getTsId(): ts.clone(ignoreAttrs=[]) for ts in tsSet if ts.getTsId() in presentTsIds}
        self.ctfDict = {ctf.getTsId(): ctf.clone(ignoreAttrs=[]) for ctf in ctfSet if ctf.getTsId() in presentTsIds}
        self.tomoDict = {tomo.getTsId(): tomo.clone() for tomo in coords.getPrecedents() if
                         tomo.getTsId() in presentTsIds}

    def convertInputStep(self):
        outPath = self._getExtraPath()
        writer = convert50_tomo.Writer()
        # Generate the particles star file
        writer.coords2Star(self.inputCoords.get(), self.tomoDict, outPath, self.coordScale.get())
        # Generate each tilt-series star file
        writer.tsSet2Star(self.tsDict, self.ctfDict, outPath)
        # Generate the tomograms star file
        writer.tomoSet2Star(self.tomoDict, self.tsDict, outPath)

    def extractSubtomos(self):
        Plugin.runRelionTomo(self, 'relion_tomo_subtomo_mpi', self.getExtractSubtomosCmd())

    def getExtractSubtomosCmd(self):
        cmd = [
            f'--p {self._getExtraPath(IN_PARTICLES_STAR)}',
            f'--t {self._getExtraPath(IN_TOMOS_STAR)}',
            f'--o {self._getExtraPath()}',
            f'--b {self.boxSize.get()}',
            f'--crop {self.croppedBoxSize.get()}',
            f'--bin {self.binningFactor.get()}',
            f'--max_dose {self.maxDose.get()}',
            f'--min_frames {self.minNoFrames.get()}'
        ]
        if self.write2dStacks.get():
            cmd.append('--stack2d')
        if self.outputInFloat16.get():
            cmd.append('--float16')
        # TODO: the native command contains a --j 10, but there's no param on its form to add a value to it...
        return ' '.join(cmd)

    #
    # def relionImportTomograms(self):
    #     Plugin.runRelionTomo(self, 'relion_tomo_import_tomograms', self._genImportTomosCmd())
    #
    # def relionImportParticles(self):
    #     Plugin.runRelionTomo(self, 'relion_tomo_import_particles', self._genImportSubtomosCmd())
    #
    # def createOutputStep(self):
    #     # Pseudosubtomos
    #     coords=self.inputCoords.get()
    #     coordSize = coords.getBoxSize()
    #     tsSamplingRate = self.tsSet.getSamplingRate()
    #     fiducialSize = int((coordSize * coords.getSamplingRate()) / (2*10))  # Radius in nm
    #
    #     psubtomoSet = createSetOfRelionPSubtomograms(self._getPath(),
    #                                                  self._getExtraPath(OPTIMISATION_SET_STAR),
    #                                                  self.inputCoords,
    #                                                  template=PSUBTOMOS_SQLITE,
    #                                                  tsSamplingRate=tsSamplingRate,
    #                                                  relionBinning=1,  # Coords are re-sampled to fit the TS size
    #                                                  boxSize=coordSize)
    #     psubtomoSet.setCoordinates3D(self.inputCoords)
    #     # Fill the set with the generated particles
    #     readSetOfPseudoSubtomograms(psubtomoSet)
    #
    #     self._defineOutputs(**{outputObjects.relionParticles.name: psubtomoSet})
    #     self._defineSourceRelation(self.inputCoords.get(), psubtomoSet)
    #     self._defineSourceRelation(self.inputCtfTs.get(), psubtomoSet)
    #
    #     # Generate the fiducial model
    #     projections = generateProjections(self._getStarFilename(OUT_PARTICLES_STAR),
    #                                       self._getStarFilename(OUT_TOMOS_STAR))
    #
    #     fiducialModelGaps = self._createSetOfLandmarkModels(suffix='Gaps')
    #     fiducialModelGaps.copyInfo(self.tsSet)
    #     fiducialModelGaps.setSetOfTiltSeries(self.inputTS) # Use the pointer better when scheduling
    #
    #     pos = 0
    #     for ts in self.tsSet:
    #         tsId = ts.getTsId()
    #
    #         if tsId not in self.matchingTSIds:
    #             continue
    #         landmarkModelGapsFilePath = os.path.join(self._getExtraPath(),
    #                                                  str(tsId) + "_gaps.sfid")
    #
    #         landmarkModelGaps = tomoObj.LandmarkModel(tsId=tsId,
    #                                                   tiltSeriesPointer=ts,
    #                                                   fileName=landmarkModelGapsFilePath,
    #                                                   modelName=None,
    #                                                   size=fiducialSize,
    #                                                   applyTSTransformation=False)
    #         landmarkModelGaps.setTiltSeries(ts)
    #
    #         # Get the projections for the tilt series
    #         tsProjections = projections[tsId]
    #
    #         for projection in tsProjections:
    #             tiltIm = projection[1] + 1
    #             chainId = projection[2] + 1
    #             xCoor = int(round(projection[3]))
    #             yCoor = int(round(projection[4]))
    #             landmarkModelGaps.addLandmark(xCoor, yCoor, tiltIm,
    #                                           chainId, 0, 0)
    #             pos += 1
    #         fiducialModelGaps.append(landmarkModelGaps)
    #
    #     self._defineOutputs(**{outputObjects.projected2DCoordinates.name: fiducialModelGaps})
    #     self._defineSourceRelation(self.tsSet, fiducialModelGaps)
    #
    # # -------------------------- INFO functions -------------------------------
    # def _warnings(self):
    #     warnMsg = []
    #     if not (self.inputTS.get().hasAlignment() and not self.inputTS.get().interpolated()):
    #         warnMsg.append('The introduced tilt series do not have an alignment transformation associated.')
    #     return warnMsg
    #
    # def _summary(self):
    #     msg = []
    #     if self.isFinished():
    #         if self.coordScale.get():
    #             msg.append('Coordinates were scaled using an scale factor of *%.2f* to be expressed considering the '
    #                        'size of the introduced tilt series' % self.coordScale.get())
    #     return msg
    #
    # # --------------------------- UTILS functions -----------------------------
    # def _genImportTomosCmd(self):
    #     acq = self.tsSet.getAcquisition()
    #     cmd = '--i %s ' % self._getStarFilename(IN_TOMOS_STAR)
    #     cmd += '--o %s ' % self._getStarFilename(OUT_TOMOS_STAR)
    #     cmd += '--hand %s ' % self._decodeHandeness()
    #     cmd += '--angpix %s ' % self.tsSet.getSamplingRate()
    #     cmd += '--voltage %s ' % acq.getVoltage()
    #     cmd += '--Cs %s ' % acq.getSphericalAberration()
    #     cmd += '--Q0 %s ' % acq.getAmplitudeContrast()
    #     if self.flipYZ.get():
    #         cmd += '--flipYZ '
    #     if self.flipZ.get():
    #         cmd += '--flipZ '
    #
    #     return cmd
    #
    # def _genImportSubtomosCmd(self):
    #     cmd = '--i %s ' % self._getStarFilename(IN_COORDS_STAR)
    #     cmd += '--o %s ' % self._getExtraPath()
    #     cmd += '--t %s ' % self._getStarFilename(OUT_TOMOS_STAR)
    #     if self.flipZCoords.get():
    #         cmd += '--flipZ '
    #     return cmd
    #
    # def _getStarFilename(self, fName):
    #     return self._getExtraPath(fName)
    #
    # def _decodeHandeness(self):
    #     return -1 if self.handeness.get() else 1
    #
    # def _simulateETomoFiles(self, imgSet, tomoSizeDict, tomoShiftDict, whiteList=None, **kwargs):
    #     """Simulate the etomo files that serve as entry point to relion4
    #     """
    #     for ts in imgSet:
    #         tsId = ts.getTsId()
    #         # Get the size of the corresponding tomogram, as it may differ from one to another
    #         tomoIdMatchDict = tomoSizeDict.get(tsId, None)
    #         if tomoIdMatchDict and (whiteList is None or tsId in whiteList):
    #             kwargs[THICKNESS] = tomoIdMatchDict[THICKNESS]
    #             kwargs[DIMS] = (tomoIdMatchDict[X_SIZE], tomoIdMatchDict[Y_SIZE])
    #
    #             # Reconstruction shifts
    #             kwargs[SHIFTS] = tomoShiftDict.get(tsId, (0,0))
    #
    #             # creating a folder where all data will be generate
    #             folderName = self._getTmpPath(tsId)
    #             makePath(folderName)
    #             # Create a symbolic link to the tiltseries image file
    #             os.symlink(os.path.abspath(ts.getFirstItem().getFileName()),
    #                        os.path.join(folderName, ts.getTsId() + '.st'))
    #             ts.writeImodFiles(folderName, **kwargs)
