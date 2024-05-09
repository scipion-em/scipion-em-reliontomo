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
from pyworkflow.utils import Message
from reliontomo import Plugin
from reliontomo.convert.convert50_tomo import getProjMatrixList, StarFileIterator, PARTICLES_TABLE, RLN_TOMONAME, \
    RLN_CENTEREDCOORDINATEXANGST, RLN_CENTEREDCOORDINATEYANGST, RLN_CENTEREDCOORDINATEZANGST
from reliontomo.objects import createSetOfRelionPSubtomograms, RelionSetOfPseudoSubtomograms
from reliontomo.constants import (IN_TOMOS_STAR, OPTIMISATION_SET_STAR, PSUBTOMOS_SQLITE, IN_PARTICLES_STAR,
                                  OUT_PARTICLES_STAR)
from reliontomo.convert import readSetOfPseudoSubtomograms, convert50_tomo
from reliontomo.protocols.protocol_re5_base_extract_subtomos_and_rec_particle import (
    ProtRelion5ExtractSubtomoAndRecParticleBase)
from reliontomo.protocols.protocol_re5_base_import_from_star import IS_RE5_PICKING_ATTR
from tomo.objects import LandmarkModel, SetOfLandmarkModels


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
        self.coordScale = None
        self.tsDict = dict()
        self.ctfDict = dict()
        self.tomoDict = dict()
        self.isRe5Picking = None

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
        coords = self.inputCoords.get()
        tsSet = self.inputTS.get()
        ctfSet = self.inputCtfTs.get()
        self.isRe5Picking = Boolean(getattr(coords, IS_RE5_PICKING_ATTR, Boolean(False).get()))

        # The ccordinates need to be re-scaled to be at the same size of the tilt-series
        self.coordScale = Float(coords.getSamplingRate() / tsSet.getSamplingRate())

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
        coordSet = self.inputCoords.get()
        # Generate the particles star file
        writer.coords2Star(coordSet, self.tomoDict, outPath,
                           coordsScale=self.coordScale.get(),
                           isRe5Picking=self.isRe5Picking)
        # Generate each tilt-series star file
        writer.tsSet2Star(self.tsDict, self.ctfDict, outPath)
        # Generate the tomograms star file
        writer.tomoSet2Star(self.tomoDict, self.tsDict, outPath, handedness=self._decodeHandedness())

    def extractSubtomos(self):
        Plugin.runRelionTomo(self, 'relion_tomo_subtomo_mpi', self.getExtractSubtomosCmd(),
                             numberOfMpi=self.numberOfMpi.get())

    def createOutputStep(self):
        tsPointer = self.inputTS
        tsSet = tsPointer.get()
        tsSRate = tsSet.getSamplingRate()
        coordsPointer = self.inputCoords
        coords = coordsPointer.get()
        psubtomoSet = createSetOfRelionPSubtomograms(self._getPath(),
                                                     self._getExtraPath(OPTIMISATION_SET_STAR),
                                                     coordsPointer,
                                                     template=PSUBTOMOS_SQLITE,
                                                     tsSamplingRate=tsSRate,
                                                     relionBinning=self.binningFactor.get(),
                                                     boxSize=coordsPointer.get().getBoxSize(),
                                                     are2dStacks=self.write2dStacks.get())
        psubtomoSet.setCoordinates3D(coordsPointer)
        # Fill the set with the generated particles
        readSetOfPseudoSubtomograms(psubtomoSet)

        # FIDUCIALS ####################################################################################################
        fiducialSize = int((coords.getBoxSize() * coords.getSamplingRate()) / (2 * 10))  # Radius in nm
        fiducialModelGaps = SetOfLandmarkModels.create(self.getPath(),
                                                       template='setOfLandmarks%s.sqlite',
                                                       suffix='Gaps')
        fiducialModelGaps.copyInfo(tsSet)
        fiducialModelGaps.setSetOfTiltSeries(tsPointer)  # Use the pointer better when scheduling
        starData = Table()
        starData.read(self._getExtraPath(OUT_PARTICLES_STAR), tableName=PARTICLES_TABLE)
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
            tsProjectionsList = getProjMatrixList(self._getExtraPath(tsId + '.star'), tomo, ts)
            for particleRow in StarFileIterator(starData, RLN_TOMONAME, tsId):
                particleCoords = np.array(
                    [self.coordScale.get() * particleRow.get(RLN_CENTEREDCOORDINATEXANGST) / tomoSRate,
                     self.coordScale.get() * particleRow.get(RLN_CENTEREDCOORDINATEYANGST) / tomoSRate,
                     self.coordScale.get() * particleRow.get(RLN_CENTEREDCOORDINATEZANGST) / tomoSRate,
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
        self._defineSourceRelation(coordsPointer, psubtomoSet)
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
            if self.coordScale.get():
                coordsFromRelion5 = self.isRe5Picking.get()
                msg.append('Coordinates were scaled using an scale factor of *%.2f* to be expressed considering the '
                           'size of the introduced tilt series' % (1 if coordsFromRelion5 else self.coordScale.get()))
        return msg

    # --------------------------- UTILS functions -----------------------------
    def getExtractSubtomosCmd(self):
        cmd = [
            f'--p {self._getExtraPath(IN_PARTICLES_STAR)}',
            f'--t {self._getExtraPath(IN_TOMOS_STAR)}',
            f'--o {self._getExtraPath()}',
            f'--max_dose {self.maxDose.get()}',
            f'--min_frames {self.minNoFrames.get()}',
            self._genCommonExtractAndRecCmd()
        ]
        if self.write2dStacks.get():
            cmd.append('--stack2d')
        if self.outputInFloat16.get():
            cmd.append('--float16')
        cmd.append('--theme classic')
        return ' '.join(cmd)

    def _decodeHandedness(self):
        return -1 if self.handedness.get() else 1
