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
import csv
from typing import Dict, Union

from emtable import Table
from pwem import ALIGN_NONE
from pwem.convert.headers import fixVolume
from pwem.emlib.image import ImageHandler
from pwem.objects import Transform
from pyworkflow.object import Integer, Float
from pyworkflow.utils import getParentFolder, yellowStr, createLink, makePath, prettyTime
from relion.convert import OpticsGroups
from reliontomo.constants import *
import numpy as np
from os.path import join, basename
from reliontomo.convert.convertBase import (checkSubtomogramFormat,
                                            getTransformInfoFromCoordOrSubtomo,
                                            WriterTomo, ReaderTomo, getTransformMatrixFromRow)
from reliontomo.objects import RelionPSubtomogram
from tomo.constants import BOTTOM_LEFT_CORNER, TR_RELION, SCIPION
from tomo.objects import Coordinate3D, SubTomogram, TomoAcquisition, Tomogram, TiltSeries, CTFTomoSeries

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
batchruntomo.a.align.AngleOffset=%(AngleOffset)f
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
        # TODO: for now, we assume the TS and the CTF have been previously well matched.
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
            oddTi, evenTi = ti.getOddEven()
            defocusU = ctfTomo.getDefocusU()
            defocusV = ctfTomo.getDefocusV()
            trMatrix = ti.getTransform().getMatrix()
            sxAngst = trMatrix[0, 2] * sRate
            syAngst = trMatrix[1, 2] * sRate
            tsTable.addRow(
                'frames/TS_01_038_-57.0.mrc',  # 1, rlnMicrographMovieName
                1,  # 2, rlnTomoTiltMovieFrameCount
                tiltAngle,  # 3, rlnTomoNominalStageTiltAngle
                tiltAxisAngle,  # 4, rlnTomoNominalTiltAxisAngle
                acqTi.getDoseInitial(),  # 5, rlnMicrographPreExposure
                # TODO: it has to be read from the mdoc from label TargetDefocus
                -4,  # 6, rlnTomoNominalDefocus
                # TODO: manage this
                'MotionCorr/job002/frames/TS_01_038_-57_0_PS.mrc',  # 7, rlnCtfPowerSpectrum
                # TODO: the tilt-images are expexted to be unstacked, as there is no index field
                oddTi,  # 8, rlnMicrographNameEven
                evenTi,  # 9, rlnMicrographNameOdd
                ti.getFileName(),  # 10, rlnMicrographName
                # TODO: check if it's used for other calculations apart from the Bayesian polishing
                'MotionCorr/job002/frames/TS_01_038_-57_0.star',  # 11, rlnMicrographMetadata
                # TODO: check if it's used for other calculations apart from the Bayesian polishing. If True, we would
                #  have to add it to our data model
                0,  # 12, rlnAccumMotionTotal
                0,  # 13, rlnAccumMotionEarly
                0,  # 14, rlnAccumMotionLate
                # TODO: check if it's used for something or if it's just an internal metric
                'CtfFind/job003/frames/TS_01_038_-57_0_PS.ctf:mrc',  # 15, rlnCtfImage
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
                0,  # 22, rlnCtfIceRingDensity
                # TODO: I've only seen this off tilt axis estimated in EMAN...
                0,  # 23, rlnTomoXTilt
                tiltAngle,  # 24, rlnTomoYTilt
                tiltAxisAngle,  # 25, rlnTomoZRot
                sxAngst,  # 26, rlnTomoXShiftAngst
                syAngst,  # 27, rlnTomoYShiftAngst
                # TODO: do we have this?
                0.5,  # 28, rlnCtfScalefactor
                )
        # Write the STAR file
        tsTable.write(getTsStarFile(ts.getTsId(), outPath))

    @staticmethod
    def tomoSet2Star(tomoDict: Dict[str, Tomogram],
                     tsDict: Dict[str, TiltSeries],
                     outPath: str,
                     unbinnedPixSize: Union[float, None] = None) -> None:
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
                    'markerDiameter': 10,  # TODO: we need the fiducial diameter at this point. Add it to model, form advanced param?
                    'rotationAngle': acq.getTiltAxisAngle(),
                    'angleOffset': 0  # TODO: check if what it is, if it's used, and how to get it from our data.
                }
                eTomoEdf = join(outPath, tsId + '.edf')
                writeEtomoEdf(eTomoEdf, eTomoDirectiveFileDict)

                # TODO: in our case, rlnMicrographOriginalPixelSize and rlnTomoTiltSeriesPixelSize are expected to be
                tomoTable.addRow(
                    tsId,  # 1, rlnTomoName
                    acq.getVoltage(),  # 2, rlnVoltage
                    acq.getSphericalAberration(),  # 3, rlnSphericalAberration
                    acq.getAmplitudeContrast(),  # 4, rlnAmplitudeContrast
                    unbinnedSRate,  # 5, rlnMicrographOriginalPixelSize
                    # TODO: add handedness to our data model
                    -1,  # 6, rlnTomoHand
                    # TODO: check if this may vary (seems not to...)
                    1,  # 7, rlnOpticsGroupName
                    tsSRate,  # 8, rlnTomoTiltSeriesPixelSize
                    getTsStarFile(tsId, outPath),  # 9, rlnTomoTiltSeriesStarFile
                    eTomoEdf,  # 10, rlnEtomoDirectiveFile
                    tomo.getSamplingRate(),  # 11, rlnTomoTomogramBinning
                    tomoX * tomoScaleFactor,  # 12, rlnTomoSizeX
                    tomoY * tomoScaleFactor,  # 13, rlnTomoSizeY
                    tomoZ * tomoScaleFactor,  # 14, rlnTomoSizeZ
                    tomoFName,  # 15, rlnTomoReconstructedTomogram
                )
        # Write the STAR file
        starFile = join(outPath, OUT_TOMOS_STAR)
        tomoTable.write(starFile)

    def coordinates2Star(self, coordSet, subtomosStar, whitelist, sRate=1, coordsScale=1):
        """Input coordsScale is used to scale the coordinates, so they are expressed in bin 1, as expected by Relion 4"""
        tomoTable = Table(columns=self._getCoordinatesStarFileLabels())
        i = 0
        for coord in coordSet.iterCoordinates():
            tsId = coord.getTomoId()

            if tsId not in whitelist:
                continue

            angles, shifts = getTransformInfoFromCoordOrSubtomo(coord, coordSet.getSamplingRate())
            # Add row to the table which will be used to generate the STAR file
            tomoTable.addRow(
                tsId,  # 1 _rlnTomoName
                coord.getObjId(),  # 2 _rlnTomoParticleId
                coord.getGroupId() if coord.getGroupId() else 1,  # 3 _rlnTomoManifoldIndex
                # coord in pix at scale of bin1
                coord.getX(RELION_3D_COORD_ORIGIN) * coordsScale,  # 4 _rlnCoordinateX
                coord.getY(RELION_3D_COORD_ORIGIN) * coordsScale,  # 5 _rlnCoordinateY
                coord.getZ(RELION_3D_COORD_ORIGIN) * coordsScale,  # 6 _rlnCoordinateZ
                # pix * Å/pix = [shifts in Å]
                shifts[0],  #* sRate,  # 7 _rlnOriginXAngst
                shifts[1],  #* sRate,  # 8 _rlnOriginYAngst
                shifts[2],  #* sRate,  # 9 _rlnOriginZAngst
                # Angles in degrees
                angles[0],  # 10 _rlnAngleRot
                angles[1],  # 11 _rlnAngleTilt
                angles[2],  # 12 _rlnAnglePsi
                # Extended fields
                int(getattr(coord, '_classNumber', -1)),  # 13_rlnClassNumber
                # Alternated 1 and 2 values
                int(getattr(coord, '_randomSubset', (i % 2) + 1)),  # 14 _rlnRandomSubset
                coord.getX(SCIPION),  # 15 _sciXCoord
                coord.getY(SCIPION),  # 16 _sciYCoord
                coord.getZ(SCIPION),  # 17 _sciZCoord
                coord.getGroupId()  # 18 _sciGroupId
            )
            i += 1
        # Write the STAR file
        tomoTable.write(subtomosStar)

    def pseudoSubtomograms2Star(self, pSubtomoSet, outStar, withPriors=False):

        logger.info("Generating particles file (%s) from pseudosubtomogram set." % outStar)
        sRate = pSubtomoSet.getSamplingRate()
        hasCoords = pSubtomoSet.getFirstItem().hasCoordinate3D()
        tomoTable = Table(columns=self._getPseudoSubtomogramStarFileLabels(hasCoords, withPriors=withPriors))

        # Write the STAR file
        optGroup = OpticsGroups.fromString(pSubtomoSet.getAcquisition().opticsGroupInfo.get())
        with open(outStar, 'w') as f:
            optGroup.toStar(f)
            # Write header first
            partsWriter = Table.Writer(f)
            partsWriter.writeTableName(PARTICLES_TABLE)
            partsWriter.writeHeader(tomoTable.getColumns())
            for pSubtomo in pSubtomoSet.iterSubtomos():
                angles, shifts = getTransformInfoFromCoordOrSubtomo(pSubtomo, pSubtomo.getSamplingRate())
                pSubtomoFile = pSubtomo.getFileName()
                pSubtomoFile = pSubtomoFile.replace(':' + MRC, '') if pSubtomoFile else FILE_NOT_FOUND
                pSubtomoCtfFile = pSubtomo.getCtfFile() if pSubtomo.getCtfFile() else FILE_NOT_FOUND

                # Add row to the table which will be used to generate the STAR file
                rowsValues = [
                    pSubtomo.getTsId(),  # _rlnTomoName #1
                    pSubtomo.getObjId(),  # rlnTomoParticleId #2
                    pSubtomo.getManifoldIndex(),  # rlnTomoManifoldIndex #3
                    # Coords in pixels
                    pSubtomo.getX(),  # _rlnCoordinateX #4
                    pSubtomo.getY(),  # _rlnCoordinateY #5
                    pSubtomo.getZ(),  # _rlnCoordinateZ #6
                    # pix * Å/pix = [shifts in Å]
                    shifts[0],  #* sRate,  # _rlnOriginXAngst #7
                    shifts[1],  #* sRate,  # _rlnOriginYAngst #8
                    shifts[2],  #* sRate,  # _rlnOriginZAngst #9
                    # Angles in degrees
                    angles[0],  # _rlnAngleRot #10
                    angles[1],  # _rlnAngleTilt #11
                    angles[2],  # _rlnAnglePsi #12

                    pSubtomo.getClassId(),  # _rlnClassNumber #13
                    pSubtomo.getRdnSubset(),  # _rlnRandomSubset #14
                ]

                if withPriors:
                    rowsValues += [angles[0], angles[1], angles[2]]
                if hasCoords:
                    rowsValues += [pSubtomo.getCoordinate3D().getX(SCIPION),  # _sciXCoord #15
                                   pSubtomo.getCoordinate3D().getY(SCIPION),  # _sciYCoord #16
                                   pSubtomo.getCoordinate3D().getZ(SCIPION),  # _sciZCoord #17
                                   pSubtomo.getCoordinate3D().getGroupId(),  # _sciGroupId #18
                                   ]

                rowsValues += [pSubtomo.getRe4ParticleName(),  # _rlnTomoParticleName #18
                               pSubtomo.getOpticsGroupId(),  # _rlnOpticsGroup #19
                               pSubtomoFile,  # _rlnImageName #20
                               pSubtomoCtfFile]  # _rlnCtfImage #21

                partsWriter.writeRowValues(rowsValues)

    def subtomograms2Star(self, subtomoSet, subtomosStar):
        logger.info("Writing relion4 star file (%s) from subtomograms." % subtomosStar)
        tomoTable = Table(columns=self.starHeaders)
        sRate = subtomoSet.getSamplingRate()
        extraPath = join(getParentFolder(subtomosStar), 'extra')
        for subtomo in subtomoSet.iterSubtomos():
            checkSubtomogramFormat(subtomo, extraPath)
            angles, shifts = getTransformInfoFromCoordOrSubtomo(subtomo, subtomo.getSamplingRate())
            ctfFile = getattr(subtomo, '_ctfImage', None)
            if ctfFile:
                ctfFile = ctfFile.get()
            classNumber = subtomo.getClassId()
            rlnClassNumber = classNumber if classNumber else 1
            rlnTomoName = subtomo.getVolName()
            rlnImageName = subtomo.getFileName().replace(':' + MRC, '')
            rlnCtfImage = ctfFile if ctfFile else FILE_NOT_FOUND
            # Coords in pixels
            rlnCoordinateX = 0
            rlnCoordinateY = 0
            rlnCoordinateZ = 0

            if subtomo.hasCoordinate3D():
                rlnCoordinateX = subtomo.getCoordinate3D().getX(BOTTOM_LEFT_CORNER)
                rlnCoordinateY = subtomo.getCoordinate3D().getY(BOTTOM_LEFT_CORNER)
                rlnCoordinateZ = subtomo.getCoordinate3D().getZ(BOTTOM_LEFT_CORNER)

            rlnAngleRot = angles[0]
            rlnAngleTilt = angles[1]
            rlnAnglePsi = angles[2]
            # pix * Å/pix = [shifts in Å]
            rlnOriginX = shifts[0]  #* sRate
            rlnOriginY = shifts[1]  #* sRate
            rlnOriginZ = shifts[2]  #* sRate
            # Angles in degrees
            rlnTiltPrior = subtomo._tiltPriorAngle.get() if hasattr(subtomo, '_tiltPriorAngle') else rlnAngleTilt
            rlnPsiPrior = subtomo._psiPriorAngle.get() if hasattr(subtomo, '_psiPriorAngle') else rlnAnglePsi
            # Add row to the table which will be used to generate the STAR file
            fieldsToAdd = [rlnTomoName,
                           rlnImageName,
                           rlnCtfImage,
                           rlnCoordinateX,
                           rlnCoordinateY,
                           rlnCoordinateZ,
                           rlnOriginX,
                           rlnOriginY,
                           rlnOriginZ,
                           rlnAngleRot,
                           rlnAngleTilt,
                           rlnTiltPrior,
                           rlnAnglePsi,
                           rlnPsiPrior,
                           rlnClassNumber]

            tomoTable.addRow(*fieldsToAdd)

        # Write the STAR file
        tomoTable.write(subtomosStar)

    @staticmethod
    def _getTomogramStarFileLabels():
        return [
            TOMO_NAME,
            TILT_SERIES_NAME,
            CTFPLOTTER_FILE,
            IMOD_DIR,
            FRACTIONAL_DOSE,
            ACQ_ORDER_FILE,
            CULLED_FILE
        ]

    @staticmethod
    def _getCoordinatesStarFileLabels(hasCoords=True, withPriors=False):
        starFileLabels = [
            TOMO_NAME,
            TOMO_PARTICLE_ID,
            MANIFOLD_INDEX,
            COORD_X,
            COORD_Y,
            COORD_Z,
            SHIFTX_ANGST,
            SHIFTY_ANGST,
            SHIFTZ_ANGST,
            ROT,
            TILT,
            PSI,
            CLASS_NUMBER,
            RANDOM_SUBSET

        ]

        if withPriors:
            starFileLabels += [ROT_PRIOR, TILT_PRIOR, PSI_PRIOR]

        if hasCoords:
            starFileLabels += [SCIPION_COORD_X, SCIPION_COORD_Y, SCIPION_COORD_Z, SCIPION_COORD_GROUP_ID]

        return starFileLabels

    @staticmethod
    def _getPseudoSubtomogramStarFileLabels(hasCoords=True, withPriors=False):
        pSubtomosLabels = Writer._getCoordinatesStarFileLabels(hasCoords, withPriors=withPriors)
        pSubtomosLabels.extend([
            TOMO_PARTICLE_NAME,
            OPTICS_GROUP,
            SUBTOMO_NAME,
            CTF_IMAGE
        ])
        return pSubtomosLabels

    @staticmethod
    def _getCtfPlotterFile(tsId, ctfPlotterDir):
        return join(ctfPlotterDir, tsId, tsId + '.defocus')

    @staticmethod
    def _genOrderListFile(prot, ts):
        """The order file expected by Relion is A 2-column, comma-separated file with the frame-order list
        of the tilt series, where the first column is the frame (image) number (starting at 1) and the second
        column is the tilt angle (in degrees).
        :param prot: current protocol object
        :param ts: TiltSeries object"""
        outputFilename = prot._getExtraPath(ts.getTsId() + '_order_list.csv')
        tiList = [ti.clone() for ti in ts]
        ind = np.argsort([ti.getAcquisitionOrder() for ti in tiList])  # Indices to get the data sorted by acqOrder
        with open(outputFilename, mode='w') as acqOrderFile:
            acqOrderFileWriter = csv.writer(acqOrderFile, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
            acqOrderList = [ti.getAcquisitionOrder() for ti in tiList]
            if min(acqOrderList) == 0:
                [acqOrderFileWriter.writerow([tiList[i].getAcquisitionOrder() + 1, tiList[i].getTiltAngle()]) for i in
                 ind]
            else:
                [acqOrderFileWriter.writerow([tiList[i].getAcquisitionOrder(), tiList[i].getTiltAngle()]) for i in ind]

        return outputFilename

    @staticmethod
    def _genCulledFileName(prot, tsId):
        return prot._getExtraPath(tsId + '_culled.mrc:mrc')


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
            x = row.get(COORD_X, 0)
            y = row.get(COORD_Y, 0)
            z = row.get(COORD_Z, 0)
            coordinate3d.setVolume(vol)
            coordinate3d.setX(float(x) * factor, RELION_3D_COORD_ORIGIN)
            coordinate3d.setY(float(y) * factor, RELION_3D_COORD_ORIGIN)
            coordinate3d.setZ(float(z) * factor, RELION_3D_COORD_ORIGIN)
            coordinate3d.setMatrix(getTransformMatrixFromRow(row, sRate=sRate), convention=TR_RELION)
            coordinate3d.setGroupId(row.get(MANIFOLD_INDEX, 1))
            # Extended fields
            coordinate3d._classNumber = Integer(row.get(CLASS_NUMBER, -1))
            coordinate3d._randomSubset = Integer(row.get(RANDOM_SUBSET, 1))

        return coordinate3d, tomoId

    def starFile2Coords3D(self, coordsSet, precedentsSet, scaleFactor):
        precedentIdDict = {}
        for tomo in precedentsSet:
            precedentIdDict[tomo.getTsId()] = tomo.clone()

        nonMatchingTomoIds = ''
        for row in self.dataTable:
            # Consider that there can be coordinates in the star file that does not belong to any of the tomograms
            # introduced
            coord, tomoId = self.gen3dCoordFromStarRow(row, coordsSet.getSamplingRate(),
                                                       precedentIdDict, factor=scaleFactor)
            if coord:
                coordsSet.append(coord)
            else:
                if tomoId not in nonMatchingTomoIds:
                    nonMatchingTomoIds += '%s ' % tomoId

        if nonMatchingTomoIds:
            logger.info(yellowStr('The star file contains coordinates that belong to tomograms not present '
                                  'in the introduced set of tomograms: %s' % nonMatchingTomoIds))

    def starFile2PseudoSubtomograms(self, outputSet):
        sRate = outputSet.getSamplingRate()
        listOfFilesToFixVolume = []
        for counter, row in enumerate(self.dataTable):
            t = Transform()
            particleFile = row.get(SUBTOMO_NAME, None)
            ctfFile = row.get(CTF_IMAGE, None)
            psubtomo = RelionPSubtomogram(fileName=particleFile,
                                          samplingRate=sRate,
                                          ctfFile=ctfFile,
                                          tsId=row.get(TOMO_NAME),
                                          classId=row.get(CLASS_NUMBER, -1),
                                          x=row.get(COORD_X),
                                          y=row.get(COORD_Y),
                                          z=row.get(COORD_Z),
                                          rdnSubset=row.get(RANDOM_SUBSET, counter % 2),  # 1 and 2 alt. by default
                                          re4ParticleName=row.get(TOMO_PARTICLE_NAME),
                                          opticsGroupId=row.get(OPTICS_GROUP, 1),
                                          manifoldIndex=row.get(MANIFOLD_INDEX, 1 if counter % 2 else -1),  # 1 and -1
                                          logLikeliCont=row.get(LOG_LIKELI_CONTRIB, -1),
                                          maxValProbDist=row.get(MAX_VALUE_PROB_DISTRIB, -1),
                                          noSignifSamples=row.get(NO_SIGNIFICANT_SAMPLES, -1)
                                          )

            # Keeping particle id
            psubtomo.setObjId(row.get(TOMO_PARTICLE_ID))

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

            # Add the files to the list of files whose header has to be corrected to be interpreted as volumes
            if particleFile is not None and particleFile != FILE_NOT_FOUND:
                listOfFilesToFixVolume.append(particleFile)
            if ctfFile is not None and ctfFile != FILE_NOT_FOUND:
                listOfFilesToFixVolume.append(ctfFile)
            # Add current pseudosubtomogram to the output set
            outputSet.append(psubtomo)

        # Keep the number of particles to compare sizes in case of subset
        outputSet.setNReParticles(self.dataTable.size())
        # Fix volume headers
        if listOfFilesToFixVolume:
            fixVolume(listOfFilesToFixVolume)

    def starFile2SubtomogramsImport(self, subtomoSet, coordSet, linkedSubtomosDir, starFilePath):
        samplingRate = subtomoSet.getSamplingRate()
        precedentsDict = {tomo.getTsId(): tomo for tomo in coordSet.getPrecedents()}
        for row, coordinate3d in zip(self.dataTable, coordSet):
            subtomo = SubTomogram()
            transform = Transform()
            origin = Transform()

            # Files
            tomoId = row.get(TOMO_NAME, FILE_NOT_FOUND)
            tomoName = precedentsDict[tomoId].getFileName()
            subtomoName = row.get(SUBTOMO_NAME, FILE_NOT_FOUND)
            tomoFolder = join(linkedSubtomosDir, tomoId)
            makePath(tomoFolder)
            linkedSubtomoName = join(tomoFolder, basename(subtomoName))
            createLink(subtomoName, linkedSubtomoName)  # Link the subtomos to the extra folder

            # Subtomograms
            tiltPrior = row.get(TILT_PRIOR, 0)
            psiPrior = row.get(PSI_PRIOR, 0)
            transform.setMatrix(coordinate3d.getMatrix(convention=TR_RELION))

            subtomo.setVolName(tomoName)
            subtomo.setFileName(linkedSubtomoName)
            subtomo.setCoordinate3D(coordinate3d)
            subtomo.setVolId(coordinate3d.getVolId())
            subtomo.setTransform(transform)
            subtomo.setAcquisition(TomoAcquisition())
            subtomo.setClassId(row.get(CLASS_NUMBER, 0))
            subtomo.setSamplingRate(samplingRate)
            subtomo._tiltPriorAngle = Float(tiltPrior)
            subtomo._psiPriorAngle = Float(psiPrior)
            subtomo.setOrigin(origin)

            # Add current subtomogram to the output set
            subtomoSet.append(subtomo)

        # Set the set of coordinates which corresponds to the current set of subtomograms
        subtomoSet.setCoordinates3D(coordSet)
