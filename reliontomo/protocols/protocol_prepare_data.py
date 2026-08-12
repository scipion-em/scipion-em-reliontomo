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
from pyworkflow.object import Float
from pyworkflow.protocol import PointerParam, BooleanParam, LEVEL_ADVANCED
from pyworkflow.utils import makePath, Message
from reliontomo import Plugin
from reliontomo.objects import createSetOfRelionPSubtomograms, RelionSetOfPseudoSubtomograms
from reliontomo.constants import (IN_TOMOS_STAR, OUT_TOMOS_STAR, IN_COORDS_STAR,
                                  OPTIMISATION_SET_STAR, OUT_PARTICLES_STAR, PSUBTOMOS_SQLITE)
from reliontomo.convert import writeSetOfTomograms, writeSetOfCoordinates, readSetOfPseudoSubtomograms
from reliontomo.protocols.protocol_base_relion import IS_RELION_50
from reliontomo.utils import generateProjections
import tomo.objects as tomoObj
from tomo.protocols.protocol_base import ProtTomoBase

# Other constants
DEFOCUS = 'defocus'
THICKNESS = 'thickness'
X_SIZE = 'x'
Y_SIZE = 'y'
DIMS = 'dims'
SHIFTS = 'shift'


class outputObjects(Enum):
    relionParticles = RelionSetOfPseudoSubtomograms
    projected2DCoordinates = tomoObj.SetOfLandmarkModels


class ProtRelionPrepareData(EMProtocol, ProtTomoBase):
    """
    Prepares tomography data and particle coordinate information for
    compatibility with Relion 4 subtomogram analysis workflows. The
    protocol converts aligned tilt-series data, tomogram geometry,
    coordinate information, and CTF metadata into the organizational
    structure expected by Relion tomography tools.

    AI Generated:

    Prepare Data for Relion 4 (ProtRelionPrepareData) — User Manual
        Overview

        The Prepare Data for Relion 4 protocol serves as a bridge
        between tomographic reconstruction workflows and Relion-based
        subtomogram processing. Its primary goal is to transform
        coordinates, tilt-series metadata, and tomographic acquisition
        information into a format that can be directly interpreted by
        Relion tomography pipelines.

        In practical cryo-electron tomography workflows, particle
        coordinates are often generated after tomogram reconstruction
        and visualization using external software environments such as
        IMOD or Scipion. However, Relion requires a precise definition
        of tilt geometry, coordinate orientation, reconstruction
        handedness, and imaging parameters before particles can be
        imported and reconstructed as pseudo-subtomograms. This
        protocol centralizes these preparation steps and ensures that
        the biological particle positions remain correctly associated
        with their originating tilt series.

        Inputs and Biological Context

        The protocol requires three major categories of information:
        three-dimensional particle coordinates, aligned tilt series,
        and CTF estimation data. Together, these datasets define the
        geometric and optical context necessary for subtomogram
        extraction and refinement.

        The coordinate set represents the biological targets of
        interest, such as ribosomes, membrane complexes, viral spikes,
        or macromolecular assemblies identified within reconstructed
        tomograms. These coordinates are linked to their parent
        tomograms and must correspond correctly to the associated tilt
        series. Any mismatch between coordinates and tilt-series
        identifiers prevents meaningful reconstruction and should be
        corrected before processing.

        The tilt series should represent the original aligned images
        used during tomogram reconstruction. Non-interpolated aligned
        tilt series are generally preferred because they preserve the
        original sampling characteristics and avoid interpolation
        artifacts that may affect downstream subtomogram refinement.

        CTF tomo series provide the defocus information required for
        Relion tomography processing. Accurate CTF metadata is
        especially important for high-resolution subtomogram averaging,
        where small inaccuracies in optical modeling may strongly
        influence refinement quality.

        Coordinate Scaling and Sampling Consistency

        One of the most important biological considerations in this
        protocol is the relationship between tomogram sampling and
        tilt-series sampling. Coordinates generated on reconstructed
        tomograms may not directly correspond to the original pixel
        size of the tilt images. The protocol automatically reconciles
        these differences so that particle locations are expressed in
        the coordinate system expected by Relion.

        This rescaling process is essential when tomograms were
        reconstructed using binning or reduced sampling rates. Without
        proper scaling, extracted subtomograms may appear shifted,
        incorrectly centered, or entirely misplaced relative to the
        biological structure of interest.

        Coordinate System Orientation

        Cryo-electron tomography workflows frequently involve changes
        in axis orientation during reconstruction and visualization.
        Different software packages may rotate, flip, or reorder axes
        in distinct ways. This protocol provides several options to
        reconcile these coordinate-system differences before Relion
        import.

        The flip YZ option is commonly required because many tomogram
        reconstruction workflows rotate the tomogram volume after
        reconstruction for visualization purposes. If this geometric
        transformation is not reproduced consistently during particle
        import, the resulting subtomograms may appear mirrored or
        incorrectly oriented.

        Similarly, flipping the Z axis may be necessary depending on
        the reconstruction convention used by the acquisition pipeline.
        Incorrect Z orientation is one of the most common causes of
        apparent particle inversion or misplaced subtomograms during
        Relion processing.

        The optional swapping of X and Y dimensions is mainly intended
        for compatibility troubleshooting in difficult reconstruction
        environments. Biological users typically leave this option
        unchanged unless particle orientations or projection geometry
        appear inconsistent after import.

        Handedness and Tilt Geometry

        Correct handedness is fundamental in cryo-electron tomography.
        The protocol allows the user to specify whether focus decreases
        or increases along the Z direction. This setting defines the
        handedness convention used during Relion import and directly
        affects how the three-dimensional geometry is interpreted.

        Incorrect handedness can produce mirrored structures,
        inconsistent particle orientations, or biologically impossible
        reconstructions. In practice, users should carefully verify the
        handedness convention of their acquisition and reconstruction
        pipeline before large-scale processing.

        Generation of Relion-Compatible Metadata

        The protocol prepares all metadata required for Relion
        tomography import procedures, including tomogram geometry,
        reconstruction dimensions, optical parameters, and coordinate
        associations. This preparation allows Relion to interpret the
        experimental geometry consistently across all tilt series.

        In addition, the protocol generates pseudo-subtomogram
        particle sets suitable for subsequent subtomogram refinement,
        classification, and averaging workflows. These pseudo-
        subtomograms preserve the biological coordinate relationships
        while adapting them to Relion conventions.

        Fiducial Projection Models

        The protocol also creates projected landmark models derived
        from the imported particle positions. These projected
        coordinates can be useful for validating geometric consistency
        across tilt images and for visual inspection of particle
        trajectories within the tomographic acquisition geometry.

        From a biological perspective, these projection models provide
        an intuitive method for verifying whether particles are
        correctly associated with the expected structures throughout
        the tilt series. Misaligned projections may indicate problems
        with coordinate scaling, tilt geometry, or tomogram alignment.

        Practical Recommendations

        Before execution, users should verify that the coordinate set,
        tilt series, and CTF information all refer to the same
        experimental dataset and share consistent identifiers. The
        protocol only processes tilt series that are common across all
        required inputs.

        In most biological workflows, the default axis orientation
        options are appropriate when tomograms were reconstructed and
        visualized using standard IMOD conventions. However, if the
        resulting subtomograms appear inverted, mirrored, or displaced,
        users should revisit the coordinate transformation options and
        test alternative configurations.

        Accurate sampling rates are particularly important. Even small
        inconsistencies between tomogram and tilt-series pixel sizes
        may produce significant coordinate shifts during subtomogram
        extraction.

        Final Perspective

        Preparing tomography data for Relion subtomogram analysis is
        not merely a file-conversion task but a critical geometric and
        biological validation step. Correct handling of coordinate
        scaling, handedness, reconstruction orientation, and imaging
        metadata directly determines whether downstream subtomogram
        refinement will produce biologically meaningful structures.

        Careful verification of coordinate consistency and tilt-series
        geometry before large-scale refinement substantially improves
        the reliability of subtomogram averaging and reduces the risk
        of difficult-to-diagnose orientation errors later in the
        workflow.
    """
    _label = 'Prepare data for Relion 4'
    _possibleOutputs = outputObjects

    def __init__(self, **args):
        EMProtocol.__init__(self, **args)
        self.tsSet = None
        self.tomoSet = None
        self.tsReducedSet = None  # Reduced list of tiltseries with only those in the coordinates
        self.matchingTSIds = None  # Unique tomogram identifiers volId used in the coordinate set. In case is a subset
        self.coordScale = Float(1)

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
        form.addParam('handeness', BooleanParam,
                      label='Does focus decrease with Z distance?',
                      default=True,
                      expertLevel=LEVEL_ADVANCED,
                      help='It is the handedness of the tilt geometry and it is used to describe '
                           'whether the focus increases or decreases as a function of Z distance.')
        form.addParam('inputTS', PointerParam,
                      help="Tilt series with alignment (non interpolated) used in the tomograms reconstruction.",
                      pointerClass='SetOfTiltSeries',
                      label="Tilt series",
                      important=True,
                      allowsNull=False)
        form.addParam('flipZCoords', BooleanParam,
                      label='Flip Z coordinate?',
                      default=False,
                      help='This option is generally False if your coordinates are displayed correctly in Scipion. '
                           'You may want to check this to True only if you see that the extracted subtomograms'
                           ' are wrong.',
                      expertLevel=LEVEL_ADVANCED)
        form.addParam('flipYZ', BooleanParam,
                      label='Has tomogram been flipped along Y and Z?',
                      default=True,
                      help='This option is generally True if the slices of your tomogram are displayed on slice Z in '
                           'Imod. '
                           'Usually, a tomogram is flipped along Y and Z (i.e. rotated around X with 90 degrees) '
                           'after the reconstruction and before the particles have been picked. This '
                           'will tell Relion to apply the same transformation to the coordinate system.')
        form.addParam('flipZ', BooleanParam,
                      label='Has the Z axis been flipped?',
                      default=True,
                      help='This option is generally True when you apply reconstrucion in Imod. This is usually used '
                           'together with the flipYZ option.')
        form.addParam('swapXY', BooleanParam,
                      label='Swap X with Y dimensions of the tilt series',
                      default=False,
                      expertLevel=LEVEL_ADVANCED,
                      help='This may be a trial and error parameter. Depending of the reconstruction path of '
                           'your tomograms we you may need to deactivate this option to get good results. '
                           'This option will be deprecated in the future')

    # -------------------------- INSERT steps functions -----------------------
    def _insertAllSteps(self):
        self._initialize()
        self._insertFunctionStep(self.convertInputStep, needsGPU=False)
        self._insertFunctionStep(self.relionImportTomograms, needsGPU=False)
        self._insertFunctionStep(self.relionImportParticles, needsGPU=False)
        self._insertFunctionStep(self.createOutputStep, needsGPU=False)

    # -------------------------- STEPS functions ------------------------------
    def _initialize(self):
        defocusDir = self._getExtraPath(DEFOCUS)
        if not exists(defocusDir):  # It can exist in case of Continue execution
            mkdir(defocusDir)
        self.coords = self.inputCoords.get()
        self.tsSet = self._getTiltSeriesNonInterpolated()

        # Compute matching TS id among coordinates and tilt series, both could be a subset
        coordsTsIds = self.coords.getTSIds()
        self.info("TsIds present in coordinates are: %s" % coordsTsIds)

        tsIds = self.tsSet.getTSIds()
        self.info("TsIds present in Tilt series are: %s" % tsIds)

        # Intersection
        self.matchingTSIds = set(coordsTsIds).intersection(tsIds)
        if len(self.matchingTSIds):
            self.info("Tilt series associated in both, coordinates and tilt series are %s" % self.matchingTSIds)
        else:
            raise Exception("There isn't any common tilt series among the coordinates and tilt series chosen.")

        self.tomoSet = self.coords.getPrecedents()
        self.inputCtfs = self.inputCtfTs.get()
        # If coordinates are referred to a set of tomograms, they'll be rescaled
        # to be expressed in bin 1, as the ts images
        if self.tomoSet:
            self.coordScale.set(self.tomoSet.getSamplingRate() / self.tsSet.getSamplingRate())

    def _getTiltSeriesNonInterpolated(self) ->tomoObj.SetOfTiltSeries:
        return self.inputTS.get() # if self.inputTS.get() is not None else \
            #getNonInterpolatedTsFromRelations(self.inputCoords.get(), self)

    def convertInputStep(self):
        # Generate defocus files
        for ctfTomo in self.inputCtfs:
            tsId = ctfTomo.getTsId()
            if tsId not in self.matchingTSIds:
                continue
            defocusPath = self._getExtraPath(DEFOCUS, ctfTomo.getTsId())
            if not exists(defocusPath):  # It can exist in case of mode Continue execution
                mkdir(defocusPath)
            generateDefocusIMODFileFromObject(ctfTomo,
                                              join(defocusPath, ctfTomo.getTsId() + '.' + DEFOCUS),
                                              isRelion=True)
        # Thickness of the tomogram and shifts
        tomoSizeDict = {}
        tomoShiftsDict ={}
        tomoList = [tomo.clone() for tomo in self.coords.getPrecedents()]
        coordScale = self.coordScale.get()
        for tomo in tomoList:
            tsId = tomo.getTsId()
            if tsId not in self.matchingTSIds:
                continue

            x, y, thickness = tomo.getDim()
            tomoSizeDict[tsId] = {X_SIZE: x * coordScale,
                                            Y_SIZE: y * coordScale,
                                            THICKNESS: thickness * coordScale}

            # Based on this: https://github.com/3dem/relion/blob/ver4.0/src/jaz/tomography/programs/convert_projections.cpp#L368
            # First element is X, second Z!
            shiftsAngs = tomo.getShiftsFromOrigin()

            # shifts are stored in Angstrom, we convert to tomo SR and then add the half and the to TS pixel size
            shiftX = int(((shiftsAngs[0] / tomo.getSamplingRate()) + x/2) * self.coordScale.get())
            shiftZ = int(((shiftsAngs[2] / tomo.getSamplingRate()) + thickness/2) * self.coordScale.get())
            tomoShiftsDict[tsId] = (shiftX, shiftZ)

            self.info("Shifts detected for %s are: %s" % (tomo.getTsId(), (shiftX, shiftZ)))

        # Simulate the etomo files that serve as entry point to relion4
        self._simulateETomoFiles(self.tsSet, tomoSizeDict, tomoShiftsDict, binned=1, binByFactor=self.coordScale,
                                 whiteList=self.matchingTSIds, swapDims=self.swapXY.get(), tltIgnoresExcluded=True)

        # Write the tomograms star file
        writeSetOfTomograms(self.tsSet,
                            self._getStarFilename(IN_TOMOS_STAR),
                            prot=self,
                            ctfPlotterParentDir=self._getExtraPath(DEFOCUS),
                            eTomoParentDir=self._getTmpPath(),
                            whiteList=self.matchingTSIds)
        # Write the particles star file
        writeSetOfCoordinates(self.inputCoords.get(),
                              self._getStarFilename(IN_COORDS_STAR),
                              self.matchingTSIds,
                              sRate=self.tsSet.getSamplingRate(),
                              coordsScale=self.coordScale.get())

    def relionImportTomograms(self):
        Plugin.runRelionTomo(self, 'relion_tomo_import_tomograms', self._genImportTomosCmd())

    def relionImportParticles(self):
        Plugin.runRelionTomo(self, 'relion_tomo_import_particles', self._genImportSubtomosCmd())

    def createOutputStep(self):
        # Pseudosubtomos
        coords=self.inputCoords.get()
        coordSize = coords.getBoxSize()
        tsSamplingRate = self.tsSet.getSamplingRate()
        fiducialSize = int((coordSize * coords.getSamplingRate()) / (2*10))  # Radius in nm

        psubtomoSet = createSetOfRelionPSubtomograms(self._getPath(),
                                                     self._getExtraPath(OPTIMISATION_SET_STAR),
                                                     self.inputCoords,
                                                     template=PSUBTOMOS_SQLITE,
                                                     tsSamplingRate=tsSamplingRate,
                                                     relionBinning=1,  # Coords are re-sampled to fit the TS size
                                                     boxSize=coordSize)
        psubtomoSet.setCoordinates3D(self.inputCoords)
        # Fill the set with the generated particles
        readSetOfPseudoSubtomograms(psubtomoSet, isRelion5=False)

        self._defineOutputs(**{outputObjects.relionParticles.name: psubtomoSet})
        self._defineSourceRelation(self.inputCoords.get(), psubtomoSet)
        self._defineSourceRelation(self.inputCtfTs.get(), psubtomoSet)

        # Generate the fiducial model
        projections = generateProjections(self._getStarFilename(OUT_PARTICLES_STAR),
                                          self._getStarFilename(OUT_TOMOS_STAR))

        fiducialModelGaps = self._createSetOfLandmarkModels(suffix='Gaps')
        fiducialModelGaps.copyInfo(self.tsSet)
        fiducialModelGaps.setSetOfTiltSeries(self.inputTS) # Use the pointer better when scheduling

        pos = 0
        for ts in self.tsSet:
            tsId = ts.getTsId()

            if tsId not in self.matchingTSIds:
                continue
            landmarkModelGapsFilePath = os.path.join(self._getExtraPath(),
                                                     str(tsId) + "_gaps.sfid")

            landmarkModelGaps = tomoObj.LandmarkModel(tsId=tsId,
                                                      tiltSeriesPointer=ts,
                                                      fileName=landmarkModelGapsFilePath,
                                                      modelName=None,
                                                      size=fiducialSize,
                                                      applyTSTransformation=False)
            landmarkModelGaps.setTiltSeries(ts)

            # Get the projections for the tilt series
            tsProjections = projections[tsId]

            for projection in tsProjections:
                tiltIm = projection[1] + 1
                chainId = projection[2] + 1
                xCoor = int(round(projection[3]))
                yCoor = int(round(projection[4]))
                landmarkModelGaps.addLandmark(xCoor, yCoor, tiltIm,
                                              chainId, 0, 0)
                pos += 1
            fiducialModelGaps.append(landmarkModelGaps)

        self._defineOutputs(**{outputObjects.projected2DCoordinates.name: fiducialModelGaps})
        self._defineSourceRelation(self.tsSet, fiducialModelGaps)

    # -------------------------- INFO functions -------------------------------
    @classmethod
    def isDisabled(cls):
        """ Return True if this Protocol is disabled.
        Disabled protocols will not be offered in the available protocols."""
        return True if IS_RELION_50 else False

    def _warnings(self):
        warnMsg = []
        if not (self.inputTS.get().hasAlignment() and not self.inputTS.get().interpolated()):
            warnMsg.append('The introduced tilt series do not have an alignment transformation associated.')
        return warnMsg

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

    def _simulateETomoFiles(self, imgSet, tomoSizeDict, tomoShiftDict, whiteList=None, **kwargs):
        """Simulate the etomo files that serve as entry point to relion4
        """
        for ts in imgSet:
            tsId = ts.getTsId()
            # Get the size of the corresponding tomogram, as it may differ from one to another
            tomoIdMatchDict = tomoSizeDict.get(tsId, None)
            if tomoIdMatchDict and (whiteList is None or tsId in whiteList):
                kwargs[THICKNESS] = tomoIdMatchDict[THICKNESS]
                kwargs[DIMS] = (tomoIdMatchDict[X_SIZE], tomoIdMatchDict[Y_SIZE])

                # Reconstruction shifts
                kwargs[SHIFTS] = tomoShiftDict.get(tsId, (0,0))

                # creating a folder where all data will be generate
                folderName = self._getTmpPath(tsId)
                makePath(folderName)
                # Create a symbolic link to the tiltseries image file
                os.symlink(os.path.abspath(ts.getFirstItem().getFileName()),
                           os.path.join(folderName, ts.getTsId() + '.st'))
                ts.writeImodFiles(folderName, **kwargs)
