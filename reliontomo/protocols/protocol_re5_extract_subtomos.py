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
import numpy as np
from emtable import Table
from pyworkflow.object import Boolean, Float
from pyworkflow.protocol import PointerParam, BooleanParam, LEVEL_ADVANCED, IntParam
from pyworkflow.utils import Message, createLink
from reliontomo import Plugin
from reliontomo.convert.convert50_tomo import getProjMatrixList, StarFileIterator, PARTICLES_TABLE, RLN_TOMONAME, \
    RLN_CENTEREDCOORDINATEXANGST, RLN_CENTEREDCOORDINATEYANGST, RLN_CENTEREDCOORDINATEZANGST
from reliontomo.objects import createSetOfRelionPSubtomograms, RelionSetOfPseudoSubtomograms
from reliontomo.constants import (OPTIMISATION_SET_STAR, PSUBTOMOS_SQLITE,
                                  OUT_PARTICLES_STAR, IN_TOMOS_STAR, GLOBAL_TABLE, RLN_TOMOTILT_SERIES_STAR_FILE)
from reliontomo.convert import readSetOfPseudoSubtomograms, convert50_tomo
from reliontomo.protocols.protocol_re5_base_extract_subtomos_and_rec_particle import (
    ProtRelion5ExtractSubtomoAndRecParticleBase)
from reliontomo.protocols.protocol_base_import_from_star import IS_RE5_PICKING_ATTR
from reliontomo.utils import getProgram
from tomo.objects import LandmarkModel, SetOfLandmarkModels, SetOfCoordinates3D


class outputObjects(Enum):
    relionParticles = RelionSetOfPseudoSubtomograms
    projected2DCoordinates = SetOfLandmarkModels


class ProtRelion5ExtractSubtomos(ProtRelion5ExtractSubtomoAndRecParticleBase):
    """extracts the relevant cropped areas of the tilt series images for each individual particle and saves them as
    CTF-premultiplied extracted 2D image stacks (or as 3D volumes).
    """
    _label = 'Extract subtomos'

    _possibleOutputs = outputObjects

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.inCoords = None
        self.coordsScaleFactor = None
        self.tsDict = dict()
        self.ctfDict = dict()
        self.tomoDict = dict()
        self.isRe5Picking = None
        self.isInSetOf3dCoords = None

    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        form.addSection(label=Message.LABEL_INPUT)
        form.addParam('inReParticles', PointerParam,
                      pointerClass='SetOfCoordinates3D, RelionSetOfPseudoSubtomograms',
                      label="Coordinates or Pseudo-subtomograms",
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
        form.addParam('handedness', BooleanParam,
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
        self._defineExtraParams(form)
        form.addParallelSection(threads=1, mpi=3)

    # -------------------------- INSERT steps functions -----------------------
    def _insertAllSteps(self):
        self._initialize()
        self._insertFunctionStep(self.convertInputStep)
        self._insertFunctionStep(self.extractSubtomos)
        self._insertFunctionStep(self.createOutputStep)

    # -------------------------- STEPS functions ------------------------------
    def _initialize(self):
        inParticles = self.getInputParticles()
        self.isInSetOf3dCoords = self.isInputSetOf3dCoords()
        coords = inParticles if self.isInSetOf3dCoords else inParticles.getCoordinates3D()
        tsSet = self.inputTS.get()
        ctfSet = self.inputCtfTs.get()
        self.isRe5Picking = Boolean(getattr(coords, IS_RE5_PICKING_ATTR, Boolean(False).get()))

        # The coordinates need to be re-scaled to be at the same size of the tilt-series
        self.coordsScaleFactor = Float(coords.getSamplingRate() / tsSet.getSamplingRate())

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
        self.inCoords = coords
        self.tsDict = {ts.getTsId(): ts.clone(ignoreAttrs=[]) for ts in tsSet if ts.getTsId() in presentTsIds}
        self.ctfDict = {ctf.getTsId(): ctf.clone(ignoreAttrs=[]) for ctf in ctfSet if ctf.getTsId() in presentTsIds}
        self.tomoDict = {tomo.getTsId(): tomo.clone() for tomo in coords.getPrecedents() if
                         tomo.getTsId() in presentTsIds}

    def convertInputStep(self):
        # Generate required star files
        coords = self.getInputParticles()
        if self.isInSetOf3dCoords:
            outPath = self._getExtraPath()
            writer = convert50_tomo.Writer()
            # Particles.star
            writer.coords2Star(coords, self.tomoDict, outPath,
                               coordsScale=self.coordsScaleFactor.get(),
                               isRe5Picking=self.isRe5Picking)
            # Tomograms.star
            writer.tomoSet2Star(self.tomoDict, self.tsDict, outPath, handedness=self._decodeHandedness())
            # Each tilt-series star file
            writer.tsSet2Star(self.tsDict, self.ctfDict, outPath)
        else:
            # Re-extraction: particles.star
            self.genInStarFile(are2dParticles=coords.are2dStacks())
            # In case of re-extraction, the tomograms file will exist and be stored as an attribute of the set, having
            # been updated if a new one is generated, like in the protocol bayesian polishing
            createLink(coords.getTomogramsStar(), self._getExtraPath(IN_TOMOS_STAR))
            # In case of re-extraction, the tilt-series star files will exist and their corresponding path will be
            # provided by the file tomograms.star

    def extractSubtomos(self):
        nMpi = self.numberOfMpi.get()
        Plugin.runRelionTomo(self,
                             getProgram('relion_tomo_subtomo', nMpi=nMpi),
                             self.getExtractSubtomosCmd(),
                             numberOfMpi=nMpi)

    def createOutputStep(self):
        tsPointer = self.inputTS
        tsSet = tsPointer.get()
        tsSRate = tsSet.getSamplingRate()
        psubtomoSet = createSetOfRelionPSubtomograms(self._getPath(),
                                                     self._getExtraPath(OPTIMISATION_SET_STAR),
                                                     self.inCoords,
                                                     template=PSUBTOMOS_SQLITE,
                                                     tsSamplingRate=tsSRate,
                                                     relionBinning=self.binningFactor.get(),
                                                     boxSize=self.inCoords.getBoxSize(),
                                                     are2dStacks=self.write2dStacks.get(),
                                                     acquisition=tsSet.getAcquisition())
        # Fill the set with the generated particles
        readSetOfPseudoSubtomograms(psubtomoSet)

        # FIDUCIALS ####################################################################################################
        fiducialSize = int((self.inCoords.getBoxSize() * self.inCoords.getSamplingRate()) / (2 * 10))  # Radius in nm
        fiducialModelGaps = SetOfLandmarkModels.create(self.getPath(),
                                                       template='setOfLandmarks%s.sqlite',
                                                       suffix='Gaps')
        fiducialModelGaps.copyInfo(tsSet)
        fiducialModelGaps.setSetOfTiltSeries(tsPointer)  # Use the pointer better when scheduling
        starData = Table()
        starData.read(self._getExtraPath(OUT_PARTICLES_STAR), tableName=PARTICLES_TABLE)
        if not self.isInSetOf3dCoords:
            tsStarFileDict = self.getTsStarFilesFromoTomgramsStar()
        particleCounter = 1

        for tsId, ts in self.tsDict.items():
            tomo = self.tomoDict[tsId]
            tomoSRate = tomo.getSamplingRate()
            landmarkModelGapsFilePath = self._getExtraPath(tsId + "_gaps.sfid")
            landmarkModelGaps = LandmarkModel(tsId=tsId,
                                              tiltSeriesPointer=ts,
                                              fileName=landmarkModelGapsFilePath,
                                              modelName=None,
                                              size=fiducialSize,
                                              applyTSTransformation=False)
            landmarkModelGaps.setTiltSeries(ts)
            tsStarFile = self._getExtraPath(tsId + '.star') if self.isInSetOf3dCoords else tsStarFileDict[tsId]
            tsProjectionsList = getProjMatrixList(tsStarFile, tomo, ts)
            for particleRow in StarFileIterator(starData, RLN_TOMONAME, tsId):
                particleCoords = np.array(
                    [self.coordsScaleFactor.get() * particleRow.get(RLN_CENTEREDCOORDINATEXANGST) / tomoSRate,
                     self.coordsScaleFactor.get() * particleRow.get(RLN_CENTEREDCOORDINATEYANGST) / tomoSRate,
                     self.coordsScaleFactor.get() * particleRow.get(RLN_CENTEREDCOORDINATEZANGST) / tomoSRate,
                     1])
                for tiltId, tomoProjection in enumerate(tsProjectionsList):
                    proj = tomoProjection.dot(particleCoords)
                    landmarkModelGaps.addLandmark(proj[0], proj[1], tiltId, particleCounter, 0, 0)
                particleCounter += 1

            fiducialModelGaps.append(landmarkModelGaps)

        ################################################################################################################

        # Define the outputs and the relations
        self._defineOutputs(**{outputObjects.relionParticles.name: psubtomoSet,
                               outputObjects.projected2DCoordinates.name: fiducialModelGaps})
        self._defineSourceRelation(self.inReParticles, psubtomoSet)
        self._defineSourceRelation(self.inputCtfTs, psubtomoSet)
        self._defineSourceRelation(tsPointer, psubtomoSet)
        self._defineSourceRelation(tsPointer, fiducialModelGaps)

    # -------------------------- INFO functions -------------------------------
    def _warnings(self):
        warnMsg = []
        if not (self.inputTS.get().hasAlignment() and not self.inputTS.get().interpolated()):
            warnMsg.append('The introduced tilt series do not have an alignment transformation associated.')
        return warnMsg

    def _summary(self):
        msg = []
        if self.isFinished():
            if self.isInSetOf3dCoords:
                inputStr = '*coordinates*'
                coordsFromRelion5 = self.isRe5Picking.get()
                scaleFactor = 1 if coordsFromRelion5 else self.coordsScaleFactor.get()
            else:
                inputStr = '*particles*'
                scaleFactor = (self.inReParticles.get().getRelionBinning() / self.binningFactor.get())

            msg.append('The %s or particles introduced were scaled using an scale factor of *%.2f* to be '
                       'expressed considering the size of the introduced tilt series' % (inputStr, scaleFactor))

        return msg

    # --------------------------- UTILS functions -----------------------------
    def getExtractSubtomosCmd(self):
        cmd = [
            self._genCommonExtractAndRecCmd(),
            f'--max_dose {self.maxDose.get()}',
            f'--min_frames {self.minNoFrames.get()}'
        ]
        if self.write2dStacks.get():
            cmd.append('--stack2d')
        if self.outputInFloat16.get():
            cmd.append('--float16')
        cmd.append('--theme classic')
        return ' '.join(cmd)

    def _decodeHandedness(self):
        return -1 if self.handedness.get() else 1

    def getTsStarFilesFromoTomgramsStar(self):
        """In re-extraction case, the tilt-series must be read from the tomograms.star file, as they may have been
        updated re-generated in some cases, such as in the CTF refinement or in the bayesian polishing.
        :return : a dictionary of structure {tsId: tsStarFile}.
        """
        tsStarDict = dict()
        tomoDataTable = Table()
        tomoDataTable.read(self._getExtraPath(IN_TOMOS_STAR), tableName=GLOBAL_TABLE)
        for row in tomoDataTable:
            tsId = row.get(RLN_TOMONAME)
            tsFile = row.get(RLN_TOMOTILT_SERIES_STAR_FILE)
            tsStarDict[tsId] = tsFile
        return tsStarDict
