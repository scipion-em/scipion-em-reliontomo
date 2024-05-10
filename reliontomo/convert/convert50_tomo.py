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
import logging
from typing import Dict, Union, List
from emtable import Table
from pwem import ALIGN_NONE
from pwem.emlib.image import ImageHandler
from pwem.objects import Transform
from pyworkflow.object import Float
from pyworkflow.utils import yellowStr, prettyTime
from relion.convert import OpticsGroups
from reliontomo.constants import *
import numpy as np
from os.path import join
from reliontomo.convert.convertBase import (getTransformInfoFromCoordOrSubtomo,
                                            WriterTomo, ReaderTomo, getTransformMatrixFromRow, genTransformMatrix)
from reliontomo.objects import RelionPSubtomogram, RelionSetOfPseudoSubtomograms
from tomo.constants import TR_RELION, SCIPION
from tomo.objects import Coordinate3D, Tomogram, TiltSeries, CTFTomoSeries, SetOfCoordinates3D, SetOfTomograms

logger = logging.getLogger(__name__)

# TILT_SERIES METADATA #################################################################################################

RLN_MICROGRAPH_MOVIENAME = 'rlnMicrographMovieName'
RLN_TOMO_TILT_MOVIE_FRAME_COUNT = 'rlnTomoTiltMovieFrameCount'
RLN_TOMO_NOMINAL_STAGE_TILT_ANGLE = 'rlnTomoNominalStageTiltAngle'
RLN_TOMO_NOMINAL_TILT_AXIS_ANGLE = 'rlnTomoNominalTiltAxisAngle'
RLN_MICROGRAPH_PRE_EXPOSURE = 'rlnMicrographPreExposure'
RLN_TOMO_NOMINAL_DEFOCUS = 'rlnTomoNominalDefocus'
RLN_CTF_POWER_SPECTRUM = 'rlnCtfPowerSpectrum'
RLN_MICROGRAPH_NAME_EVEN = 'rlnMicrographNameEven'
RLN_MICROGRAPH_NAME_ODD = 'rlnMicrographNameOdd'
RLN_MICROGRAPH_NAME = 'rlnMicrographName'
RLN_MICROGRAPH_METADATA = 'rlnMicrographMetadata'
RLN_ACCUM_MOTION_TOTAL = 'rlnAccumMotionTotal'
RLN_ACCUM_MOTION_EARLY = 'rlnAccumMotionEarly'
RLN_ACCUM_MOTION_LATE = 'rlnAccumMotionLate'
RLN_CTF_IMAGE = 'rlnCtfImage'
RLN_DEFOCUS_U = 'rlnDefocusU'
RLN_DEFOCUS_V = 'rlnDefocusV'
RLN_CTF_ASTIGMATISM = 'rlnCtfAstigmatism'
RLN_DEFOCUS_ANGLE = 'rlnDefocusAngle'
RLN_CTF_FIGURE_OF_MERIT = 'rlnCtfFigureOfMerit'
RLN_CTF_MAX_RESOLUTION = 'rlnCtfMaxResolution'
RLN_CTF_ICE_RING_DENSITY = 'rlnCtfIceRingDensity'
RLN_TOMO_X_TILT = 'rlnTomoXTilt'
RLN_TOMO_Y_TILT = 'rlnTomoYTilt'
RLN_TOMO_Z_ROT = 'rlnTomoZRot'
RLN_TOMO_X_SHIFT_ANGST = 'rlnTomoXShiftAngst'
RLN_TOMO_Y_SHIFT_ANGST = 'rlnTomoYShiftAngst'
RLN_CTF_SCALEFACTOR = 'rlnCtfScalefactor'

tsStarFields = [
    RLN_MICROGRAPH_MOVIENAME,
    RLN_TOMO_TILT_MOVIE_FRAME_COUNT,
    RLN_TOMO_NOMINAL_STAGE_TILT_ANGLE,
    RLN_TOMO_NOMINAL_TILT_AXIS_ANGLE,
    RLN_MICROGRAPH_PRE_EXPOSURE,
    RLN_TOMO_NOMINAL_DEFOCUS,
    RLN_CTF_POWER_SPECTRUM,
    RLN_MICROGRAPH_NAME_EVEN,
    RLN_MICROGRAPH_NAME_ODD,
    RLN_MICROGRAPH_NAME,
    RLN_MICROGRAPH_METADATA,
    RLN_ACCUM_MOTION_TOTAL,
    RLN_ACCUM_MOTION_EARLY,
    RLN_ACCUM_MOTION_LATE,
    RLN_CTF_IMAGE,
    RLN_DEFOCUS_U,
    RLN_DEFOCUS_V,
    RLN_CTF_ASTIGMATISM,
    RLN_DEFOCUS_ANGLE,
    RLN_CTF_FIGURE_OF_MERIT,
    RLN_CTF_MAX_RESOLUTION,
    RLN_CTF_ICE_RING_DENSITY,
    RLN_TOMO_X_TILT,
    RLN_TOMO_Y_TILT,
    RLN_TOMO_Z_ROT,
    RLN_TOMO_X_SHIFT_ANGST,
    RLN_TOMO_Y_SHIFT_ANGST,
    RLN_CTF_SCALEFACTOR
]

# TOMOGRAMS METADATA ###################################################################################################
GLOBAL_TABLE = 'global'

RLN_TOMONAME = 'rlnTomoName'
RLN_VOLTAGE = 'rlnVoltage'
RLN_SPHERICALABERRATION = 'rlnSphericalAberration'
RLN_AMPLITUDECONTRAST = 'rlnAmplitudeContrast'
RLN_MICROGRAPHORIGINALPIXELSIZE = 'rlnMicrographOriginalPixelSize'
RLN_TOMOHAND = 'rlnTomoHand'
RLN_OPTICSGROUPNAME = 'rlnOpticsGroupName'
RLN_TOMOTILT_SERIES_PIXEL_SIZE = 'rlnTomoTiltSeriesPixelSize'
RLN_TOMOTILT_SERIES_STAR_FILE = 'rlnTomoTiltSeriesStarFile'
RLN_ETOMO_DIRECTIVE_FILE = 'rlnEtomoDirectiveFile'
RLN_TOMOTOMOGRAM_BINNING = 'rlnTomoTomogramBinning'
RLN_TOMOSIZEX = 'rlnTomoSizeX'
RLN_TOMOSIZEY = 'rlnTomoSizeY'
RLN_TOMOSIZEZ = 'rlnTomoSizeZ'
RLN_TOMORECONSTRUCTED_TOMOGRAM = 'rlnTomoReconstructedTomogram'

tomoStarFields = [
    RLN_TOMONAME,
    RLN_VOLTAGE,
    RLN_SPHERICALABERRATION,
    RLN_AMPLITUDECONTRAST,
    RLN_MICROGRAPHORIGINALPIXELSIZE,
    RLN_TOMOHAND,
    RLN_OPTICSGROUPNAME,
    RLN_TOMOTILT_SERIES_PIXEL_SIZE,
    RLN_TOMOTILT_SERIES_STAR_FILE,
    RLN_ETOMO_DIRECTIVE_FILE,
    RLN_TOMOTOMOGRAM_BINNING,
    RLN_TOMOSIZEX,
    RLN_TOMOSIZEY,
    RLN_TOMOSIZEZ,
    RLN_TOMORECONSTRUCTED_TOMOGRAM
]

# COORDINATES METADATA #################################################################################################
PARTICLES_TABLE = 'particles'
OPTICS_TABLE = 'optics'

RLN_COORDINATEX = 'rlnCoordinateX'
RLN_COORDINATEY = 'rlnCoordinateY'
RLN_COORDINATEZ = 'rlnCoordinateZ'
RLN_TOMOSUBTOMOGRAMROT = 'rlnTomoSubtomogramRot'
RLN_TOMOSUBTOMOGRAMTILT = 'rlnTomoSubtomogramTilt'
RLN_TOMOSUBTOMOGRAMPSI = 'rlnTomoSubtomogramPsi'
RLN_ANGLEROT = 'rlnAngleRot'
RLN_ANGLETILT = 'rlnAngleTilt'
RLN_ANGLEPSI = 'rlnAnglePsi'
RLN_ANGLETILTPRIOR = 'rlnAngleTiltPrior'
RLN_ANGLEPSIPRIOR = 'rlnAnglePsiPrior'

coordsStarFields = [
    RLN_TOMONAME,
    RLN_COORDINATEX,
    RLN_COORDINATEY,
    RLN_COORDINATEZ,
    RLN_TOMOSUBTOMOGRAMROT,
    RLN_TOMOSUBTOMOGRAMTILT,
    RLN_TOMOSUBTOMOGRAMPSI,
    RLN_ANGLEROT,
    RLN_ANGLETILT,
    RLN_ANGLEPSI,
    RLN_ANGLETILTPRIOR,
    RLN_ANGLEPSIPRIOR,
    SCIPION_COORD_X,
    SCIPION_COORD_Y,
    SCIPION_COORD_Z,
    SCIPION_COORD_GROUP_ID
]

# PARTICLES METADATA ###################################################################################################

RLN_TOMOSUBTOMOGRAMROT = 'rlnTomoSubtomogramRot'
RLN_TOMOSUBTOMOGRAMTILT = 'rlnTomoSubtomogramTilt'
RLN_TOMOSUBTOMOGRAMPSI = 'rlnTomoSubtomogramPsi'
RLN_ANGLETILTPRIOR = 'rlnAngleTiltPrior'
RLN_ANGLEPSIPRIOR = 'rlnAnglePsiPrior'
RLN_OPTICSGROUP = 'rlnOpticsGroup'
RLN_TOMOPARTICLENAME = 'rlnTomoParticleName'
RLN_TOMOVISIBLEFRAMES = 'rlnTomoVisibleFrames'
RLN_CTFIMAGE = 'rlnCtfImage'
RLN_IMAGENAME = 'rlnImageName'
RLN_ORIGINXANGST = 'rlnOriginXAngst'
RLN_ORIGINYANGST = 'rlnOriginYAngst'
RLN_ORIGINZANGST = 'rlnOriginZAngst'
RLN_CENTEREDCOORDINATEXANGST = 'rlnCenteredCoordinateXAngst'
RLN_CENTEREDCOORDINATEYANGST = 'rlnCenteredCoordinateYAngst'
RLN_CENTEREDCOORDINATEZANGST = 'rlnCenteredCoordinateZAngst'
RLN_GROUPNUMBER = 'rlnGroupNumber'
RLN_CLASSNUMBER = 'rlnClassNumber'
RLN_NORMCORRECTION = 'rlnNormCorrection'
RLN_RANDOMSUBSET = 'rlnRandomSubset'
RLN_LOGLIKELYCONTRIB = 'rlnLogLikeliContribution'
RLN_MAXVALPROBDISTIRB = 'rlnMaxValueProbDistribution'
RLN_NSIGNIFSAMPLES = 'rlnNrOfSignificantSamples'

particles2dStarFields = [
    RLN_TOMONAME,
    RLN_TOMOSUBTOMOGRAMROT,
    RLN_TOMOSUBTOMOGRAMTILT,
    RLN_TOMOSUBTOMOGRAMPSI,
    RLN_ANGLEROT,
    RLN_ANGLETILT,
    RLN_ANGLEPSI,
    RLN_ANGLETILTPRIOR,
    RLN_ANGLEPSIPRIOR,
    RLN_OPTICSGROUP,
    RLN_TOMOPARTICLENAME,
    RLN_TOMOVISIBLEFRAMES,
    RLN_IMAGENAME,
    RLN_ORIGINXANGST,
    RLN_ORIGINYANGST,
    RLN_ORIGINZANGST,
    RLN_CENTEREDCOORDINATEXANGST,
    RLN_CENTEREDCOORDINATEYANGST,
    RLN_CENTEREDCOORDINATEZANGST,
    RLN_GROUPNUMBER,
    RLN_CLASSNUMBER,
    RLN_NORMCORRECTION,
    RLN_RANDOMSUBSET,
    RLN_LOGLIKELYCONTRIB,
    RLN_MAXVALPROBDISTIRB,
    RLN_NSIGNIFSAMPLES,
]

particles3dStarFields = [
    RLN_TOMONAME,
    RLN_TOMOSUBTOMOGRAMROT,
    RLN_TOMOSUBTOMOGRAMTILT,
    RLN_TOMOSUBTOMOGRAMPSI,
    RLN_ANGLEROT,
    RLN_ANGLETILT,
    RLN_ANGLEPSI,
    RLN_ANGLETILTPRIOR,
    RLN_ANGLEPSIPRIOR,
    RLN_OPTICSGROUP,
    RLN_TOMOPARTICLENAME,
    RLN_CTFIMAGE,
    RLN_IMAGENAME,
    RLN_ORIGINXANGST,
    RLN_ORIGINYANGST,
    RLN_ORIGINZANGST,
    RLN_CENTEREDCOORDINATEXANGST,
    RLN_CENTEREDCOORDINATEYANGST,
    RLN_CENTEREDCOORDINATEZANGST,
    RLN_GROUPNUMBER,
    RLN_CLASSNUMBER,
    RLN_NORMCORRECTION,
    RLN_RANDOMSUBSET,
    RLN_LOGLIKELYCONTRIB,
    RLN_MAXVALPROBDISTIRB,
    RLN_NSIGNIFSAMPLES,
]

# OTHERS
eyeMatrix3x3 = np.eye(3)

# RELION 5 3D COORDINATES
R5_ROT_ATTRIB = '_rlnAngleRot'
R5_TILT_ATTRIB = '_rlnAngleTilt'
R5_PSI_ATTRIB = '_rlnAnglePsi'
R5_TILT_PRIOR_ATTRIB = '_rlnAngleTiltPrior'
R5_PSI_PRIO_ATTRIB = '_rlnAnglePsiPrior'


def getTsStarFile(tsId: str, outPath: str) -> str:
    """
    It generates the star file name of a given tsId.
    :param tsId: tilt-series identifier.
    :param outPath: path (only path, no filename) in which the star file will be generated.
    """
    return join(outPath, tsId + '.star')


# TODO: update these lines with the ones generated by relion5 and check the values to be updated, concretely
#  batchruntomo.a.align.AngleOffset=12.91
def writeEtomoEdf(fn, paramsDict):
    template = """
#%(date)s
Setup.DataSource=CCD
ReconstructionState.InvalidEdgeFunctionsA=no result
Setup.BackupDirectory=
ReconstructionState.InvalidEdgeFunctionsB=no result
ProcessTrack.PostProcessing=Not started
Setup.Combine.ManualCleanup=false
Setup.AxisA.TiltAngle.Type=Extract
Setup.FiducialessAlignmentB=false
Setup.FiducialessAlignmentA=false
ProcessTrack.FinalAlignedStack-A=Not started
ProcessTrack.FinalAlignedStack-B=Not started
Setup.tiltalign.TargetPatchSizeXandY=700,700
Setup.AxisB.TiltAngle.RangeMin=%(minTilt)f
Setup.Combine.TempDirectory=
Setup.ImageFile.ImageFilenameStyle=MRC
Setup.Combine.UseList=
Setup.Combine.PatchBoundaryYMax=0
Setup.WholeTomogramSampleA=true
Setup.DatasetName=%(name)s
batchruntomo.MadeZFactorsA=false
batchruntomo.UsedLocalAlignmentsA=false
batchruntomo.a.align.AxisZShift=0.0
batchruntomo.a.align.AngleOffset=%(angleOffset)f
batchruntomo.a.SeedingDone=true
batchruntomo.Track.A.LightBeads=0
batchruntomo.OrigImageStackExt=mrc
Setup.FiducialDiameter=%(markerDiameter)f
Setup.WholeTomogramSampleB=true
Setup.SetFEIPixelSize=false
ReconstructionState.A.AdjustOrigin=true
ProcessTrack.RevisionNumber=2.0
Setup.Combine.PatchBoundaryYMin=0
Setup.Combine.MaxPatchBoundaryZMax=0
ProcessTrack.CoarseAlignment-A=Not started
Setup.ViewType=Single View
ProcessTrack.CoarseAlignment-B=Not started
ReconstructionState.TrimvolFlipped=no result
Setup.Track.B.TrackMethod=Seed
Setup.Track.A.TrackMethod=Seed
Setup.A.SizeToOutputInXandY=/
Setup.FinalStack.b.UseExpandCircleIterations=true
ProcessTrack.TomogramGeneration-A=Not started
ProcessTrack.TomogramGeneration-B=Not started
ProcessTrack.CleanUp=Not started
Setup.AxisA.TiltAngle.RangeMin=%(minTilt)f
Setup.AxisB.TiltAngle.TiltAngleFilename=
Setup.Pos.A.NewDialog=true
Setup.Squeezevol.LinearInterpolation=false
Setup.Stack.B.Is.Twodir=false
Setup.Combine.RevisionNumber=1.2
Setup.Stack.B.CTF.AutoFit.RangeAndStep=-Infinity,-Infinity
Setup.UseLocalAlignmentsB=true
Setup.UseLocalAlignmentsA=true
Setup.AxisB.TiltAngle.RangeStep=1.0
Setup.Combine.FiducialMatchListA=
Setup.Combine.FiducialMatchListB=
Setup.B.SizeToOutputInXandY=/
Setup.FinalStack.a.UseExpandCircleIterations=true
Setup.Combine.ModelBased=false
ReconstructionState.MadeZFactorsB=no result
Setup.Version.Etomo.Created=4.11.24
ReconstructionState.MadeZFactorsA=no result
Setup.FinalStackBinningB=8
ReconstructionState.SqueezevolFlipped=no result
Setup.FinalStackBinningA=8
ProcessTrack.FiducialModel-B=Not started
ProcessTrack.FiducialModel-A=Not started
Setup.Combine.PatchBoundaryXMax=0
ProcessTrack.Setup=Complete
ProcessTrack.PreProcessing-A=Not started
ProcessTrack.PreProcessing-B=Not started
ProcessTrack.FineAlignment-B=Not started
ProcessTrack.FineAlignment-A=Not started
Setup.AxisA.TiltAngle.TiltAngleFilename=
Setup.OrigRawImageStackExt=mrc
Setup.AxisA.TiltAngle.RangeStep=1.0
Setup=-Infinity,-Infinity
Setup.Combine.PatchBoundaryXMin=0
Setup.RevisionNumber=1.12
Setup.Track.B.SeedModel.Transfer=true
Setup.Track.A.Raptor.UseRawStack=false
Setup.PixelSize=%(pixelSize)f
ReconstructionState.B.AdjustOrigin=true
ReconstructionState.Combine.ScriptsCreated=no result
Setup.Combine.Transfer=true
ReconstructionState.UsedLocalAlignmentsA=no result
ReconstructionState.UsedLocalAlignmentsB=no result
Setup.Combine.FiducialMatch=BothSides
Setup.AxisType=Single Axis
Setup.ImageRotationA=%(rotationAngle)f
ProcessTrack.TomogramPositioning-A=Not started
Setup.ImageRotationB=
Setup.Stack.A.Is.Twodir=false
Setup.Pos.B.NewDialog=true
ProcessTrack.TomogramPositioning-B=Not started
Setup.Combine.PatchBoundaryZMax=0
Setup.DefaultGpuProcessing=false
Setup.Track.A.SeedModel.Auto=true
Setup.RawImageStackExt=mrc
Setup.Combine.PatchSize=M
Setup.AxisB.TiltAngle.Type=Extract
Setup.Combine.PatchBoundaryZMin=0
Setup.FinalStack.b.ExpandCircleIterations=3
Setup.DefaultParallel=false
Setup.Version.Etomo.Modified=4.11.24
Setup.tiltalign.NumberOfLocalPatchesXandY=5,5
Setup.Combine.PatchRegionModel=
ReconstructionState.NewstFiducialessAlignmentA=no result
ReconstructionState.NewstFiducialessAlignmentB=no result
Setup.FinalStack.a.ExpandCircleIterations=3
ProcessTrack.TomogramCombination=Not started
        """
    with open(fn, 'w') as f:
        f.write(template % paramsDict)


class Writer(WriterTomo):
    """ Helper class to convert from Scipion SetOfTomograms and SetOfSubTomograms to star files ."""

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    def tsSet2Star(self,
                   tsDict: Dict[str, TiltSeries],
                   ctfDict: Dict[str, CTFTomoSeries],
                   outPath: str) -> None:
        """
        It generates a tilt-series star file for each tilt-series contained in the given tsDict that matches via tsId
        with a CTFTomoSeries from ctfDict.
        :param tsDict: dictionary of type {tsId: TiltSeries}
        :param ctfDict: dictionary of type {tsId: CTFTomoSeries}
        :param outPath: path (only path, no filename) in which the star file will be generated.
        """
        for tsId, ts in tsDict.items():
            ctf = ctfDict.get(tsId, None)
            if ctf:
                self.ts2Star(ts, ctf, outPath)

    @staticmethod
    def ts2Star(ts, ctf, outPath):
        """
        It writes a tilt-series star file in relion5 format. Fields (output of the command execution
        relion_refine --print_metadata_labels):

        rlnMicrographMovieName #1 (string) : Name of a micrograph movie stack
        rlnTomoTiltMovieFrameCount #2 (int)    : Number of frames in the tilt series movies
        rlnTomoNominalStageTiltAngle #3 (double) : Nominal value for the stage tilt angle
        rlnTomoNominalTiltAxisAngle #4 (double) : Nominal value for the angle of the tilt axis
        rlnMicrographPreExposure #5 (double) : Pre-exposure dose in electrons per square Angstrom
        rlnTomoNominalDefocus #6 (double) : Nominal value for the defocus in the tilt series image
        rlnCtfPowerSpectrum #7 (string) : Power spectrum for CTF estimation
        rlnMicrographNameEven #8 (string) : Micrograph summed from even frames of motion corrected movie
        rlnMicrographNameOdd #9 (string) : Micrograph summed from odd frames of motion corrected movie
        rlnMicrographName #10 (string) : Name of a micrograph
        rlnMicrographMetadata #11 (string) : Name of a micrograph metadata file
        rlnAccumMotionTotal #12 (double) : Accumulated global motion during the entire movie (in A)
        rlnAccumMotionEarly #13 (double) : Accumulated global motion during the first frames of the movie (in A)
        rlnAccumMotionLate #14 (double) : Accumulated global motion during the last frames of the movie (in A)
        rlnCtfImage #15 (string) : Name of an image with all CTF values
        rlnDefocusU #16 (double) : Defocus in U-direction (in Angstroms, positive values for underfocus)
        rlnDefocusV #17 (double) : Defocus in V-direction (in Angstroms, positive values for underfocus)
        rlnCtfAstigmatism #18 (double) : Absolute value of the difference between defocus in U- and V-direction (in A)
        rlnDefocusAngle #19 (double) : Angle between X and defocus U direction (in degrees)
        rlnCtfFigureOfMerit #20 (double) : Figure of merit for the fit of the CTF (not used inside relion_refine)
        rlnCtfMaxResolution #21 (double) : Estimated maximum resolution (in A) of significant CTF Thon rings
        rlnCtfIceRingDensity #22 (double) : Power of the image in the ice ring frequency range (0.25-0.28 A-1)
        rlnTomoXTilt #23 (double) : Euler angle for rotation of tomogram around X-axis
        rlnTomoYTilt #24 (double) : Euler angle for rotation of tomogram around Y-axis
        rlnTomoZRot #25 (double) : Euler angle for rotation of tomogram around Z-axis
        rlnTomoXShiftAngst #26 (double) : X-translation (in A) to align the projection of a tomogram with the tilt series image
        rlnTomoYShiftAngst #27 (double) : Y-translation (in A) to align the projection of a tomogram with the tilt series image
        rlnCtfScalefactor #28 (double) : Linear scale-factor on the CTF (values between 0 and 1)

        Example:
            frames/TS_01_038_-57.0.mrc            8    -56.99850    85.000000   114.000000     -4.00000
            MotionCorr/job002/frames/TS_01_038_-57_0_PS.mrc MotionCorr/job002/frames/TS_01_038_-57_0_EVN.mrc
            MotionCorr/job002/frames/TS_01_038_-57_0_ODD.mrc MotionCorr/job002/frames/TS_01_038_-57_0.mrc
            MotionCorr/job002/frames/TS_01_038_-57_0.star     6.904044     0.000000     6.904044
            CtfFind/job003/frames/TS_01_038_-57_0_PS.ctf:mrc 38313.808594 38228.539062    85.269531    54.530354
            -0.08450    16.288582     0.711630     0.000000    -57.00000    85.080380    37.591033   126.622629
            0.544639
        """
        # TODO: for now, we skip the exclusion stuff, leave it for the future once it works with the simplest case
        # TODO: Tilt angle and tilt axis angle are stored separately for the initial and refined values. It is logical
        #  to think that the initial ones are there for tracking from the upper part of the pipeline and that the ones
        #  used are the refined ones.
        sRate = ts.getSamplingRate()
        tiltAxisAngle = ts.getAcquisition().getTiltAxisAngle()
        tsTable = Table(columns=tsStarFields)
        for ti, ctfTomo in zip(ts, ctf):
            acqTi = ti.getAcquisition()
            tiltAngle = ti.getTiltAngle()
            oddEven = ti.getOddEven()
            if oddEven:
                oddTi, evenTi = oddEven
            else:
                oddTi = FILE_NOT_FOUND
                evenTi = FILE_NOT_FOUND

            defocusU = ctfTomo.getDefocusU()
            defocusV = ctfTomo.getDefocusV()
            trMatrix = ti.getTransform().getMatrix() if ti.getTransform() is not None else eyeMatrix3x3
            iTrMatrix = np.linalg.inv(trMatrix)
            rotAngle = np.rad2deg(np.arccos(trMatrix[0, 0]))
            sxAngst = iTrMatrix[0, 2] * sRate
            syAngst = iTrMatrix[1, 2] * sRate
            tsTable.addRow(
                FILE_NOT_FOUND,  # 1, rlnMicrographMovieName
                1,  # 2, rlnTomoTiltMovieFrameCount
                tiltAngle,  # 3, rlnTomoNominalStageTiltAngle
                tiltAxisAngle,  # 4, rlnTomoNominalTiltAxisAngle
                acqTi.getDoseInitial(),  # 5, rlnMicrographPreExposure
                # TODO: it has to be read from the mdoc from label TargetDefocus
                -1.5,  # 6, rlnTomoNominalDefocus
                # TODO: manage this
                FILE_NOT_FOUND,  # 7, rlnCtfPowerSpectrum
                # TODO: the tilt-images are expexted to be unstacked, as there is no index field
                oddTi,  # 8, rlnMicrographNameEven
                evenTi,  # 9, rlnMicrographNameOdd
                f'{ti.getIndex()}@{ti.getFileName()}:mrcs',  # 10, rlnMicrographName
                # TODO: check if it's used for other calculations apart from the Bayesian polishing
                FILE_NOT_FOUND,  # 11, rlnMicrographMetadata
                # TODO: check if it's used for other calculations apart from the Bayesian polishing. If True, we would
                #  have to add it to our data model
                0,  # 12, rlnAccumMotionTotal
                0,  # 13, rlnAccumMotionEarly
                0,  # 14, rlnAccumMotionLate
                # TODO: check if it's used for something or if it's just an internal metric
                # src/reconstructor.cpp
                # C++
                # ·
                # master
                # 			Image<RFLOAT> Ictf;
                # 			FileName fn_ctf;
                # 			if (!DF.getValue(EMDL_CTF_IMAGE, fn_ctf, p))
                # 				REPORT_ERROR("ERROR: cannot find rlnCtfImage for 3D CTF correction!");
                # 			Ictf.read(fn_ctf);
                # 			// If there is a redundant half, get rid of it
                ctfTomo.getPsdFile() if ctfTomo.getPsdFile() else FILE_NOT_FOUND,  # 15, rlnCtfImage
                defocusU,  # 16, rlnDefocusU
                defocusV,  # 17, rlnDefocusV
                abs(defocusU - defocusV),  # 18, rlnCtfAstigmatism
                ctfTomo.getDefocusAngle(),  # 19, rlnDefocusAngle
                # TODO: from Relion label definition: "not used inside relion_refine", but may be _fitQuality from
                #  our model
                0,  # 20, rlnCtfFigureOfMerit
                # TODO: check if CTFFind's Estimated maximum resolution (in A) of significant CTF Thon rings is the
                #  same in other plugins that estimate the CTF
                ctfTomo.getResolution(),  # 21, rlnCtfMaxResolution
                # TODO: check if this value is used and, in that case, if we have this somewhere or have to store it
                #  in the data model
                ctfTomo.getFitQuality() if ctfTomo.getFitQuality() is not None else 0,  # 22, rlnCtfIceRingDensity
                # TODO: I've only seen this off tilt axis estimated in EMAN...
                0,  # 23, rlnTomoXTilt
                tiltAngle,  # 24, rlnTomoYTilt
                rotAngle,  # 25, rlnTomoZRot
                sxAngst,  # 26, rlnTomoXShiftAngst
                syAngst,  # 27, rlnTomoYShiftAngst
                np.cos(np.deg2rad(tiltAngle)),  # 28, rlnCtfScalefactor
            )
        # Write the STAR file
        tsId = ts.getTsId()
        tsTable.write(getTsStarFile(tsId, outPath), tableName=tsId)

    @staticmethod
    def tomoSet2Star(tomoDict: Dict[str, Tomogram],
                     tsDict: Dict[str, TiltSeries],
                     outPath: str,
                     unbinnedPixSize: Union[float, None] = None,
                     handedness: int = -1) -> None:
        """
        It writes a tomograms star file in relion5 format. Only the tomograms that matches with a tilt-series from
        tsDic will be added. Fields (output of the command execution relion_refine --print_metadata_labels):

        rlnTomoName #1 (string) : Arbitrary name for a tomogram
        rlnVoltage #2 (double) : Voltage of the microscope (in kV)
        rlnSphericalAberration #3 (double) : Spherical aberration (in millimeters)
        rlnAmplitudeContrast #4 (double) : Amplitude contrast (as a fraction, i.e. 10% = 0.1)
        rlnMicrographOriginalPixelSize #5 (double) : Pixel size of original movie before binning in Angstrom/pixel.
        rlnTomoHand #6 (double) : Handedness of a tomogram (i.e. slope of defocus over the image-space z coordinate)
        rlnOpticsGroupName #7 (string) : The name of a group of particles with identical optical properties
        rlnTomoTiltSeriesPixelSize #8 (double) : Pixel size of the original tilt series
        rlnTomoTiltSeriesStarFile #9 (string) : Tilt series starfile
        rlnEtomoDirectiveFile #10 (string) : Location of the etomo directive file (.edf) from tilt series alignment
        rlnTomoTomogramBinning #11 (double) : Binning level of a  reconstructed tomogram
        rlnTomoSizeX #12 (int)    : Width of a bin-1 tomogram in pixels
        rlnTomoSizeY #13 (int)    : Height of a bin-1 tomogram in pixels
        rlnTomoSizeZ #14 (int)    : Depth of a bin-1 tomogram in pixels
        rlnTomoReconstructedTomogram #15 (string) : File name of a reconstructed tomogram

        Example:
                 TS_01   300.000000     2.700000     0.100000     0.675000     -1.00000    optics1     1.350000
                 Tomograms/job022/tilt_series/TS_01.star AlignTiltSeries/job021/external/TS_01/TS_01.edf
                 7.407407         4000         4000         2000 Tomograms/job022/tomograms/rec_TS_01.mrc

        :param tomoDict: dictionary of type {tsId: Tomogram}
        :param tsDict: dictionary of type {tsId: TiltSeries}
        :param outPath: path (only path, no filename) in which the star file will be generated.
        :param unbinnedPixSize: original pixel size. If not provided, it will be considered the same as the tilt-series
        :param handedness: handedness of the tomograms (i.e. slope of defocus over the image-space z coordinate)
        """

        # We'll use the acquisition from the corresponding tilt-series to each tomogram as it is more reliable (in case
        # of imported tomograms it may not have some of the data required here). Also, we'll do it tomo by tomo, as
        # there may be heterogeneity in some of the parameters. Thus, each tilt-series acquisition is more reliable
        # than the one for the whole set
        tomoTable = Table(columns=tomoStarFields)
        ih = ImageHandler()
        for tsId, tomo in tomoDict.items():
            ts = tsDict.get(tsId, None)
            if ts:
                tsSRate = ts.getSamplingRate()
                unbinnedSRate = unbinnedPixSize if unbinnedPixSize else tsSRate
                acq = ts.getAcquisition()
                tomoFName = tomo.getFileName()
                tomoX, tomoY, tomoZ, _ = ih.getDimensions(tomoFName)
                tomoScaleFactor = tomo.getSamplingRate() / unbinnedSRate

                eTomoDirectiveFileDict = {
                    'date': prettyTime(),
                    'name': tsId,
                    'pixelSize': tsSRate,
                    'minTilt': acq.getAngleMin(),
                    'markerDiameter': 10,
                    # TODO: we need the fiducial diameter at this point. Add it to model, form advanced param?
                    'rotationAngle': acq.getTiltAxisAngle(),
                    'angleOffset': -0.12  # TODO: check if what it is, if it's used, and how to get it from our data.
                }
                eTomoEdf = join(outPath, tsId + '.edf')
                # writeEtomoEdf(eTomoEdf, eTomoDirectiveFileDict)

                # TODO: in our case, rlnMicrographOriginalPixelSize and rlnTomoTiltSeriesPixelSize are expected to be
                tomoTable.addRow(
                    tsId,  # 1, rlnTomoName
                    acq.getVoltage(),  # 2, rlnVoltage
                    acq.getSphericalAberration(),  # 3, rlnSphericalAberration
                    acq.getAmplitudeContrast(),  # 4, rlnAmplitudeContrast
                    unbinnedSRate,  # 5, rlnMicrographOriginalPixelSize
                    handedness,  # 6, rlnTomoHand
                    # TODO: check if this may vary (seems not to...)
                    'optics1',  # 7, rlnOpticsGroupName
                    tsSRate,  # 8, rlnTomoTiltSeriesPixelSize
                    getTsStarFile(tsId, outPath),  # 9, rlnTomoTiltSeriesStarFile
                    eTomoEdf,  # 10, rlnEtomoDirectiveFile
                    tomo.getSamplingRate() / tsSRate,  # 11, rlnTomoTomogramBinning
                    tomoX * tomoScaleFactor,  # 12, rlnTomoSizeX
                    tomoY * tomoScaleFactor,  # 13, rlnTomoSizeY
                    tomoZ * tomoScaleFactor,  # 14, rlnTomoSizeZ
                    tomoFName,  # 15, rlnTomoReconstructedTomogram
                )
        # Write the STAR file
        tomoTable.write(join(outPath, IN_TOMOS_STAR), tableName=GLOBAL_TABLE)

    @staticmethod
    def coords2Star(coordSet: SetOfCoordinates3D,
                    tomoDict: Dict[str, Tomogram],
                    outPath: str,
                    isRe5Picking: bool = False,
                    coordsScale: float = 1.0
                    ):
        """
        It writes a particles star file in relion5 format. Fields (output of the command execution
        relion_refine --print_metadata_labels):

        rlnTomoName #1 (string) : Arbitrary name for a tomogram
        rlnCoordinateX #2 (double) : X-Position of an image in a micrograph (in pixels)
        rlnCoordinateY #3 (double) : Y-Position of an image in a micrograph (in pixels)
        rlnCoordinateZ #4 (double) : Z-Position of an image in a 3D micrograph, i.e. tomogram (in pixels)
        rlnTomoSubtomogramRot #5 (double) : First Euler angle of a subtomogram (rot, in degrees)
        rlnTomoSubtomogramTilt #6 (double) : Second Euler angle of a subtomogram (tilt, in degrees)
        rlnTomoSubtomogramPsi #7 (double) : Third Euler angle of a subtomogram (psi, in degrees)
        rlnAngleRot #8 (double) : First Euler angle (rot, in degrees)
        rlnAngleTilt #9 (double) : Second Euler angle (tilt, in degrees)
        rlnAnglePsi #10 (double) : Third Euler angle (psi, in degrees)
        rlnAngleTiltPrior #11 (double) : Center of the prior (in degrees) on the second Euler angle (tilt)
        rlnAnglePsiPrior #12 (double) : Center of the prior (in degrees) on the third Euler angle (psi)

        Example:
            TS_43	1798.595537	1030.038284	1208.171910	-178.979203	87.527230	-22.439980	0.000000	90.000000
            0.000000	90.000000	0.000000

        :param coordSet: SetOfCoordinates3D
        :param tomoDict: dictionary of type {tsId: Tomogram}
        :param outPath: path (only path, no filename) in which the star file will be generated.
        :param isRe5Picking: used to indicate if the coordinates were picked using the picker provided by
        Relion5. If set to True, the coordinates won't be scaled, as they're written scaled to bin1 in the generated
        star file (the one which is being imported).
        :param coordsScale: used to scale the coordiantes to the tilt-series size. Ignored and set to 1 if
        isRe5Picking = True.
        """
        particlesTable = Table(columns=coordsStarFields)
        sRate = coordSet.getSamplingRate()
        coordsScale = 1 if isRe5Picking else coordsScale
        for tsId, tomo in tomoDict.items():
            for coord in coordSet.iterCoordinates(volume=tomo):
                angles, _ = getTransformInfoFromCoordOrSubtomo(coord, sRate)
                particlesTable.addRow(
                    tsId,  # 1, rlnTomoName
                    coord.getX(RELION_3D_COORD_ORIGIN) * coordsScale,  # 2, rlnCoordinateX
                    coord.getY(RELION_3D_COORD_ORIGIN) * coordsScale,  # 3, rlnCoordinateY
                    coord.getZ(RELION_3D_COORD_ORIGIN) * coordsScale,  # 4, rlnCoordinateZ
                    angles[0],  # 5, rlnTomoSubtomogramRot
                    angles[1],  # 6, rlnTomoSubtomogramTilt
                    angles[2],  # 7, rlnTomoSubtomogramPsi
                    getattr(coord, R5_ROT_ATTRIB, Float(0)).get(),  # 8, rlnAngleRot
                    getattr(coord, R5_TILT_ATTRIB, Float(0)).get(),  # 9, rlnAngleTilt
                    getattr(coord, R5_PSI_ATTRIB, Float(0)).get(),  # 10, rlnAnglePsi
                    getattr(coord, R5_TILT_PRIOR_ATTRIB, Float(0)).get(),  # 11, rlnAngleTiltPrior
                    getattr(coord, R5_PSI_PRIO_ATTRIB, Float(0)).get(),  # 12, rlnAnglePsiPrior
                    # Scipion fields
                    coord.getX(SCIPION),  # _sciXCoord
                    coord.getY(SCIPION),  # _sciYCoord
                    coord.getZ(SCIPION),  # _sciZCoord
                    coord.getGroupId()  # _sciGroupId
                )
        # Write the STAR file
        particlesTable.write(join(outPath, IN_PARTICLES_STAR), tableName=PARTICLES_TABLE)

    @staticmethod
    def pseudoSubtomograms2Star(pSubtomoSet: RelionSetOfPseudoSubtomograms,
                                outPath: str,
                                are2dParticles: bool = True):

        logger.info("Generating particles file from pseudosubtomogram set.")
        particlesStarFields = particles2dStarFields if are2dParticles else particles3dStarFields
        sRate = pSubtomoSet.getSamplingRate()
        hasCoords = pSubtomoSet.getFirstItem().hasCoordinate3D()
        if hasCoords:
            sciCoordsFields = [
                SCIPION_COORD_X,
                SCIPION_COORD_Y,
                SCIPION_COORD_Z,
                SCIPION_COORD_GROUP_ID
            ]
            particlesStarFields.extend(sciCoordsFields)
        tomoTable = Table(columns=particlesStarFields)

        # Write the STAR file
        optGroup = OpticsGroups.fromString(pSubtomoSet.getAcquisition().opticsGroupInfo.get())
        with open(join(outPath, IN_PARTICLES_STAR), 'w') as f:
            optGroup.toStar(f)
            # Write header first
            partsWriter = Table.Writer(f)
            partsWriter.writeTableName(PARTICLES_TABLE)
            partsWriter.writeHeader(tomoTable.getColumns())
            for pSubtomo in pSubtomoSet.iterSubtomos():
                angles, shifts = getTransformInfoFromCoordOrSubtomo(pSubtomo, sRate)
                # Fields 12 and 13 is different depending on if the particles are 2D or 3D
                imageName = pSubtomo.getFileName()
                if are2dParticles:
                    field12 = pSubtomo.getVisibleFrames()
                    field13 = imageName
                else:
                    field12 = imageName
                    field13 = pSubtomo.getCtfFile() if pSubtomo.getCtfFile() else FILE_NOT_FOUND

                # Add row to the table which will be used to generate the STAR file
                rowsValues = [
                    pSubtomo.getTsId(),  # 1, rlnTomoName
                    angles[0],  # 2, rlnTomoSubtomogramRot
                    angles[1],  # 3, rlnTomoSubtomogramTilt
                    angles[2],  # 4, rlnTomoSubtomogramPsi
                    pSubtomo.getRot(),  # 5, rlnAngleRot
                    pSubtomo.getTilt(),  # 6, rlnAngleTilt
                    pSubtomo.getPsi(),   # 7, rlnAnglePsi
                    pSubtomo.getTiltPrior(),  # 8, rlnAngleTiltPrior
                    pSubtomo.getPsiPrior(),  # 9, rlnAnglePsiPrior
                    pSubtomo.getOpticsGroupId(),  # 10, rlnOpticsGroup
                    pSubtomo.getRelionParticleName(),  # 11, rlnTomoParticleName
                    field12,  # 12, rlnTomoVisibleFrames (for 2D particles) or rlnImageName (for 3D particles)
                    field13,  # 13, rlnImageName (for 2D particles) or rlnCtfImage (for 3D particles)
                    # pix * Å/pix = [shifts in Å]
                    shifts[0] * sRate,  # 14, rlnOriginXAngst
                    shifts[1] * sRate,  # 15, rlnOriginYAngst
                    shifts[2] * sRate,  # 16, rlnOriginZAngst,
                    pSubtomo.getXInImg(),  # 17, rlnCenteredCoordinateXAngst
                    pSubtomo.getYInImg(),  # 18, rlnCenteredCoordinateYAngst
                    pSubtomo.getZInImg(),  # 19, rlnCenteredCoordinateZAngst
                    pSubtomo.getGroupId(),  # 20, rlnGroupNumber
                    pSubtomo.getClassId(),  # 21, rlnClassNumber
                    pSubtomo.getNormCorrection(),  # 22, rlnNormCorrection
                    pSubtomo.getRdnSubset(),  # 23, rlnRandomSubset
                    pSubtomo.getLogLikeliContribution(),  # 24, rlnLogLikeliContribution
                    pSubtomo.getMaxValueProbDistribution(),  # 25,  getMaxValueProbDistribution
                    pSubtomo.getNrOfSignificantSamples()  # 26, rlnNrOfSignificantSamples
                ]
                if hasCoords:
                    rowsValues += [pSubtomo.getCoordinate3D().getX(SCIPION),  # 27, sciXCoord
                                   pSubtomo.getCoordinate3D().getY(SCIPION),  # 28, sciYCoord
                                   pSubtomo.getCoordinate3D().getZ(SCIPION),  # 29, sciZCoord
                                   pSubtomo.getCoordinate3D().getGroupId(),  # 30, sciGroupId
                                   ]
                partsWriter.writeRowValues(rowsValues)


class Reader(ReaderTomo):
    ALIGNMENT_LABELS = [
        SHIFTX_ANGST,
        SHIFTY_ANGST,
        SHIFTZ_ANGST,
        ROT,
        TILT,
        PSI,
    ]

    def __init__(self, starFile, dataTable, **kwargs):
        super().__init__(starFile, dataTable)
        self._shifts = np.zeros(3)
        self._angles = np.zeros(3)
        self._alignType = kwargs.get('alignType', ALIGN_NONE)
        self._pixelSize = kwargs.get('pixelSize', 1.0)

    @staticmethod
    def gen3dCoordFromStarRow(row, sRate, precedentIdDict, factor=1):
        coordinate3d = None
        tomoId = row.get(TOMO_NAME)
        vol = precedentIdDict.get(tomoId, None)
        if vol:
            coordinate3d = Coordinate3D()
            x = row.get(RLN_COORDINATEX, 0)
            y = row.get(RLN_COORDINATEY, 0)
            z = row.get(RLN_COORDINATEZ, 0)
            coordinate3d.setVolume(vol)
            coordinate3d.setX(float(x) * factor, RELION_3D_COORD_ORIGIN)
            coordinate3d.setY(float(y) * factor, RELION_3D_COORD_ORIGIN)
            coordinate3d.setZ(float(z) * factor, RELION_3D_COORD_ORIGIN)
            trMatrix = getCoordsTransformMatrixFromRow(row, sRate=sRate)
            coordinate3d.setMatrix(trMatrix, convention=TR_RELION)
            # Extended fields
            setattr(coordinate3d, R5_ROT_ATTRIB, Float(row.get(RLN_ANGLEROT, 0)))
            setattr(coordinate3d, R5_TILT_ATTRIB, Float(row.get(RLN_ANGLETILT, 0)))
            setattr(coordinate3d, R5_PSI_ATTRIB, Float(row.get(RLN_ANGLEPSI, 0)))
            setattr(coordinate3d, R5_TILT_PRIOR_ATTRIB, Float(row.get(RLN_ANGLETILTPRIOR, 0)))
            setattr(coordinate3d, R5_PSI_PRIO_ATTRIB, Float(row.get(RLN_ANGLEPSIPRIOR, 0)))

        return coordinate3d, tomoId

    def starFile2Coords3D(self,
                          coordsSet: SetOfCoordinates3D,
                          tomogramsSet: SetOfTomograms,
                          scaleFactor: bool = 1):
        """ Converts the contents of a preloaded star file into Scipion SetOfCoordinates3D.
        :param coordsSet: SetOfCoordinates3D that will be filled with the contest from the loaded star file
        :param tomogramsSet: introduced SetOfTomograms.
        :param scaleFactor: used to scale the coordinates to the size of the tomograms."""
        presentTsIds = tomogramsSet.getTSIds()
        precedentIdDict = {tomo.getTsId(): tomo.clone() for tomo in tomogramsSet if tomo.getTsId() in presentTsIds}

        nonMatchingTomoIds = ''
        for row in self.dataTable:
            # Consider that there can be coordinates in the star file that does not belong to any of the tomograms
            # introduced
            coord, tomoId = self.gen3dCoordFromStarRow(row,
                                                       coordsSet.getSamplingRate(),
                                                       precedentIdDict,
                                                       factor=scaleFactor)
            if coord:
                coordsSet.append(coord)
            else:
                if tomoId not in nonMatchingTomoIds:
                    nonMatchingTomoIds += '%s ' % tomoId

        if nonMatchingTomoIds:
            logger.info(yellowStr('The star file contains coordinates that belong to tomograms not present '
                                  'in the introduced set of tomograms: %s' % nonMatchingTomoIds))

    def starFile2PseudoSubtomograms(self, outputSet):
        """Reads the data_particles table of a generated particles.star file. (output of the command execution
        relion_refine --print_metadata_labels):

        rlnTomoName #1 (string) : Arbitrary name for a tomogram
        rlnTomoSubtomogramRot #2 (double) : First Euler angle of a subtomogram (rot, in degrees)
        rlnTomoSubtomogramTilt #3 (double) : Second Euler angle of a subtomogram (tilt, in degrees)
        rlnTomoSubtomogramPsi #4 (double) : Third Euler angle of a subtomogram (psi, in degrees)
        rlnAngleRot #5 (double) : First Euler angle (rot, in degrees)
        rlnAngleTilt #6 (double) : Second Euler angle (tilt, in degrees)
        rlnAnglePsi #7 (double) : Third Euler angle (psi, in degrees)
        rlnAngleTiltPrior #8 (double) : Center of the prior (in degrees) on the second Euler angle (tilt)
        rlnAnglePsiPrior #9 (double) : Center of the prior (in degrees) on the third Euler angle (psi)
        rlnOpticsGroup #10 (int)    : Group of particles with identical optical properties
        rlnTomoParticleName #11 (string) : Name of each individual particle
        rlnTomoVisibleFrames #12 (vector<int>) : Frames fromt he tilt series that are included in the 2D stack of a pseudo-subtomogram
        rlnImageName #13 (string) : Name of an image
        rlnOriginXAngst #14 (double) : X-coordinate (in Angstrom) for the origin of rotation
        rlnOriginYAngst #15 (double) : Y-coordinate (in Angstrom) for the origin of rotation
        rlnOriginZAngst #16 (double) : Z-coordinate (in Angstrom) for the origin of rotation
        rlnCenteredCoordinateXAngst #17 (double) : X-Position of an image in a micrograph (in Angstroms, with the center being 0,0)
        rlnCenteredCoordinateYAngst #18 (double) : Y-Position of an image in a micrograph (in Angstroms, with the center being 0,0)
        rlnCenteredCoordinateZAngst #19 (double) : Z-Position of an image in a micrograph (in Angstroms, with the center being 0,0)
        rlnGroupNumber #20 (int)    : The number of a group of images
        rlnClassNumber #21 (int)    : Class number for which a particle has its highest probability
        rlnNormCorrection #22 (double) : Normalisation correction value for an image
        rlnRandomSubset #23 (int)    : Random subset to which this particle belongs
        rlnLogLikeliContribution #24 (double) : Contribution of a particle to the log-likelihood target function
        rlnMaxValueProbDistribution #25 (double) : Maximum value of the (normalised) probability function for a particle
        rlnNrOfSignificantSamples #26 (int)    : Number of orientational/class assignments (for a particle)
        with sign.probabilities in the 1st pass of adaptive oversampling.

        Example:
             TS_01   -130.15816    85.886437   130.158156    26.433646    11.507468   167.399460     0.000000
             0.000000            1   TS_01/1 Runs/003096_ProtRelion5ExtractSubtomos/extra/Subtomograms/TS_01/1_data.mrc
             Runs/003096_ProtRelion5ExtractSubtomos/extra/Subtomograms/TS_01/1_weights.mrc    -49.87757    -32.09807
             23.812184  -1390.63902  -1580.08937   660.057917            1            1     1.000000            1
             3.950164e+06     0.409595            8
        """
        sRate = outputSet.getSamplingRate()
        for counter, row in enumerate(self.dataTable):
            t = Transform()
            particleFile = row.get(RLN_IMAGENAME, None)
            psubtomo = RelionPSubtomogram(fileName=particleFile,
                                          samplingRate=sRate,
                                          tsId=row.get(RLN_TOMONAME),
                                          classId=row.get(RLN_CLASSNUMBER, -1),
                                          # Coords. in Relion they are in angstroms, while in Scipion they're in pixels
                                          x=row.get(RLN_ORIGINXANGST, 0),
                                          y=row.get(RLN_ORIGINYANGST, 0),
                                          z=row.get(RLN_ORIGINZANGST, 0),
                                          xInImg=row.get(RLN_CENTEREDCOORDINATEXANGST, 0),
                                          yInImg=row.get(RLN_CENTEREDCOORDINATEYANGST, 0),
                                          zInImg=row.get(RLN_CENTEREDCOORDINATEZANGST, 0),
                                          rdnSubset=row.get(RANDOM_SUBSET, counter % 2),  # 1 and 2 alt. by default
                                          relionParticleName=row.get(RLN_TOMOPARTICLENAME),
                                          visibleFrames=row.get(RLN_TOMOVISIBLEFRAMES, [0, 0, 0]),
                                          ctfFile=row.get(RLN_CTFIMAGE, FILE_NOT_FOUND),
                                          opticsGroupId=row.get(OPTICS_GROUP, 1),
                                          manifoldIndex=row.get(MANIFOLD_INDEX, 1 if counter % 2 else -1),  # 1 and -1
                                          logLikeliCont=row.get(LOG_LIKELI_CONTRIB, -1),
                                          maxValProbDist=row.get(MAX_VALUE_PROB_DISTRIB, -1),
                                          noSignifSamples=row.get(NO_SIGNIFICANT_SAMPLES, -1),
                                          rot=row.get(RLN_ANGLEROT, 0),
                                          tilt=row.get(RLN_ANGLETILT, 0),
                                          psi=row.get(RLN_ANGLEPSI, 0),
                                          tiltPrior=row.get(RLN_ANGLETILTPRIOR, 0),
                                          psiPrior=row.get(RLN_ANGLEPSIPRIOR, 0),
                                          groupId=row.get(RLN_GROUPNUMBER, 1),
                                          normCorrection=row.get(RLN_NORMCORRECTION, 0),
                                          )

            # TODO: decide what to do with this
            # Keeping particle id
            # psubtomo.setObjId(row.get(TOMO_PARTICLE_ID))

            # Set the coordinate3D
            if row.get(SCIPION_COORD_X) is not None:  # Assume that the coordinates exists
                sciCoord = Coordinate3D()
                sciCoord.setX(row.get(SCIPION_COORD_X), SCIPION)
                sciCoord.setY(row.get(SCIPION_COORD_Y), SCIPION)
                sciCoord.setZ(row.get(SCIPION_COORD_Z), SCIPION)
                sciCoord.setGroupId(row.get(SCIPION_COORD_GROUP_ID, 1))
                sciCoord.setTomoId(row.get(TOMO_NAME))
                psubtomo.setCoordinate3D(sciCoord)

            # Set the transformation matrix
            t.setMatrix(getTransformMatrixFromRow(row, sRate=sRate))
            psubtomo.setTransform(t)
            psubtomo.setIndex(counter)
            psubtomo.setClassId(row.get(RLN_CLASSNUMBER, 1))

            # Add current pseudosubtomogram to the output set
            outputSet.append(psubtomo)

        # Keep the number of particles to compare sizes in case of subset
        outputSet.setNReParticles(self.dataTable.size())


def getProjMatrixList(tsStarFile: str, tomogram: Tomogram, ts: TiltSeries) -> List[np.ndarray]:
    # / *From
    # Alister
    # Burt
    # *
    # *tilt_image_center = tilt_image_dimensions / 2
    # *specimen_center = tomogram_dimensions / 2
    # *
    # *  # Transformations, defined in order of application
    # *s0 = S(-specimen_center)  # put specimen center-of-rotation at the origin
    # *r0 = Rx(euler_angles['rlnTomoXTilt'])  # rotate specimen around X-axis
    # *r1 = Ry(euler_angles['rlnTomoYTilt'])  # rotate specimen around Y-axis
    # *r2 = Rz(euler_angles['rlnTomoZRot'])  # rotate specimen around Z-axis
    # *s1 = S(specimen_shifts)  # shift projected specimen in xy (camera) plane
    # *s2 = S(tilt_image_center)  # move specimen back into tilt-image coordinate system
    # *
    # *  # compose matrices
    # *transformations = s2 @ s1 @ r2 @ r1 @ r0 @ s0
    # *
    # * /
    # specimen_shifts(xshift_angst / optics.pixelSize, yshift_angst / optics.pixelSize, 0.);
    prjMatrixList = []
    dataTable = Table()
    dataTable.read(tsStarFile)
    tsSRate = ts.getSamplingRate()
    ih = ImageHandler()
    tomoXDim, tomoYDim, tomoZDim, _ = ih.getDimensions(tomogram.getFileName())
    tsXDim, tsYDim, _, _ = ih.getDimensions(ts.getFirstItem().getFileName())
    for row in dataTable:
        xRotAngle = row.get(RLN_TOMO_X_TILT, None)
        yRotAngle = row.get(RLN_TOMO_Y_TILT, None)
        zRotAngle = row.get(RLN_TOMO_Z_ROT, None)
        sxAngst = row.get(RLN_TOMO_X_SHIFT_ANGST, None)
        syAngst = row.get(RLN_TOMO_Y_SHIFT_ANGST, None)
        # Rotate specimen around X-axis
        r0 = gen3dRotXMatrix(xRotAngle)
        r1 = gen3dRotYMatrix(yRotAngle)
        r2 = gen3dRotZMatrix(zRotAngle)
        # Translations
        s0 = genTranslationMatrix(-tomoXDim / 2, -tomoYDim / 2, -tomoZDim / 2)
        s1 = genTranslationMatrix(sxAngst / tsSRate, syAngst / tsSRate, 0)
        s2 = genTranslationMatrix(tsXDim / 2, tsYDim / 2, 0)
        # Projection matrix
        # prjMatrix = s2 @ s1 @ r2 @ r1 @ r0 @ s0
        prjMatrix = s2 @ s1 @ r2 @ r1 @ r0
        # logger.info(yellowStr(prjMatrix))
        # logger.info(f's2 =\n{s2}')
        # logger.info(f's1 =\n{s1}')
        # logger.info(f'r2 =\n{r2}')
        # logger.info(f'r1 =\n{r1}')
        # logger.info(f'r0 =\n{r0}')
        # logger.info(f's0 =\n{s0}')
        prjMatrixList.append(prjMatrix)

    return prjMatrixList


def gen3dRotXMatrix(angleInDeg):
    angleInRad = np.radians(angleInDeg)
    return np.array([[1, 0, 0, 0],
                     [0, np.cos(angleInRad), -np.sin(angleInRad), 0],
                     [0, np.sin(angleInRad), np.cos(angleInRad), 0],
                     [0, 0, 0, 1]])


def gen3dRotYMatrix(angleInDeg):
    angleInRad = np.radians(angleInDeg)
    return np.array([[np.cos(angleInRad), 0, np.sin(angleInRad), 0],
                     [0, 1, 0, 0],
                     [-np.sin(angleInRad), 0, np.cos(angleInRad), 0],
                     [0, 0, 0, 1]])


def gen3dRotZMatrix(angleInDeg):
    angleInRad = np.radians(angleInDeg)
    return np.array([[np.cos(angleInRad), -np.sin(angleInRad), 0, 0],
                     [np.sin(angleInRad), np.cos(angleInRad), 0, 0],
                     [0, 0, 1, 0],
                     [0, 0, 0, 1]])


def genTranslationMatrix(sx, sy, sz):
    trMatrix = np.eye(4)
    trMatrix[0, 3] = sx
    trMatrix[1, 3] = sy
    trMatrix[2, 3] = sz
    return trMatrix


def getCoordsTransformMatrixFromRow(row, sRate=1):
    shiftx = 0
    shifty = 0
    shiftz = 0
    rot = row.get(RLN_TOMOSUBTOMOGRAMROT, 0)
    tilt = row.get(RLN_TOMOSUBTOMOGRAMTILT, 0)
    psi = row.get(RLN_TOMOSUBTOMOGRAMPSI, 0)
    return genTransformMatrix(shiftx, shifty, shiftz, rot, tilt, psi, sRate)


class StarFileIterator:
    def __init__(self, star_data, field_name, field_value):
        self.star_data = star_data
        self.field_name = field_name
        self.field_value = field_value
        self.index = 0

    def __iter__(self):
        return self

    def __next__(self):
        while self.index < len(self.star_data):
            if self.star_data[self.index].get(self.field_name) == self.field_value:
                result = self.star_data[self.index]
                self.index += 1
                return result
            self.index += 1
        raise StopIteration
