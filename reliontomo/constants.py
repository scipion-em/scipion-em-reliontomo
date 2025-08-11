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
import json
from os.path import join

from tomo.constants import BOTTOM_LEFT_CORNER

RELION = 'relion'
V4_0 = '4.0'
RELIONTOMO_HOME = 'RELIONTOMO_HOME'
RELIONTOMO_CUDA_LIB = 'RELION_CUDA_LIB'
RELIONTOMO_DEFAULT_VERSION = V4_0
RELIONTOMO_DEFAULT = RELION + '-' + RELIONTOMO_DEFAULT_VERSION
V30_VALIDATION_MSG = 'This version of Reliontomo plugin requires Relion 3.0 binaries. ' \
                     'Please install it if necessary via scipion3 installb relion-3.0'

STAR_COMPARE_ENTRY_POINT = 'compareStar'
MRC = 'mrc'
PSUBTOMOS_SQLITE = 'relionPSubtomograms%s.sqlite'
RELION_3D_COORD_ORIGIN = BOTTOM_LEFT_CORNER

# Relion 4 tomogram star file fields
TOMO_NAME = 'rlnTomoName'
TILT_SERIES_NAME = 'rlnTomoTiltSeriesName'
CTFPLOTTER_FILE = 'rlnTomoImportCtfPlotterFile'
IMOD_DIR = 'rlnTomoImportImodDir'
FRACTIONAL_DOSE = 'rlnTomoImportFractionalDose'
ACQ_ORDER_FILE = 'rlnTomoImportOrderList'
CULLED_FILE = 'rlnTomoImportCulledFile'

# Relion 3 specific star file fields
TOMO_NAME_30 = 'rlnMicrographName'
CTF_MISSING_WEDGE = 'rlnCtfImage'
MAGNIFICATION = 'rlnMagnification'
PIXEL_SIZE = 'rlnDetectorPixelSize'
SHIFTX = 'rlnOriginX'
SHIFTY = 'rlnOriginY'
SHIFTZ = 'rlnOriginZ'

# Relion 4 subtomogram star file fields
SUBTOMO_NAME = 'rlnImageName'
COORD_X = 'rlnCoordinateX'
COORD_Y = 'rlnCoordinateY'
COORD_Z = 'rlnCoordinateZ'
SHIFTX_ANGST = 'rlnOriginXAngst'
SHIFTY_ANGST = 'rlnOriginYAngst'
SHIFTZ_ANGST = 'rlnOriginZAngst'
ROT = 'rlnAngleRot'
TILT = 'rlnAngleTilt'
PSI = 'rlnAnglePsi'
ROT_PRIOR = 'rlnAngleRotPrior'
TILT_PRIOR = 'rlnAngleTiltPrior'
PSI_PRIOR = 'rlnAnglePsiPrior'
CLASS_NUMBER = 'rlnClassNumber'

# Relion pseudo-subtomograms fields
TOMO_PARTICLE_ID = 'rlnTomoParticleId'
TOMO_PARTICLE_NAME = 'rlnTomoParticleName'
MANIFOLD_INDEX = 'rlnTomoManifoldIndex'
RANDOM_SUBSET = 'rlnRandomSubset'
OPTICS_GROUP = 'rlnOpticsGroup'
CTF_IMAGE = 'rlnCtfImage'
LOG_LIKELI_CONTRIB = 'rlnLogLikeliContribution'
MAX_VALUE_PROB_DISTRIB = 'rlnMaxValueProbDistribution'
NO_SIGNIFICANT_SAMPLES = 'rlnNrOfSignificantSamples'

# Optimisation set fields
OPT_PARTICLES_STAR = 'rlnTomoParticlesFile'
OPT_TOMOS_STAR = 'rlnTomoTomogramsFile'
OPT_TRAJECTORIES_STAR = 'rlnTomoTrajectoriesFile'
OPT_MANIFOLDS_STAR = 'rlnTomoManifoldsFile'
OPT_FSC_STAR = 'rlnTomoReferenceFscFile'

# Scipion Coordinates fields
SCIPION_COORD_X = 'sciXCoord'
SCIPION_COORD_Y = 'sciYCoord'
SCIPION_COORD_Z = 'sciZCoord'
SCIPION_COORD_GROUP_ID = 'sciGroupId'

# Warp Coordinates fields
WRP_COORDINATE_X = 'wrpCoordinateX'
WRP_COORDINATE_Y = 'wrpCoordinateY'
WRP_COORDINATE_Z = 'wrpCoordinateZ'

FILE_NOT_FOUND = 'FileNotFound'

RELION_30_TOMO_LABELS = [TOMO_NAME_30,
                         COORD_X,
                         COORD_Y,
                         COORD_Z,
                         SUBTOMO_NAME,
                         CTF_MISSING_WEDGE,
                         MAGNIFICATION,
                         PIXEL_SIZE,
                         ROT,
                         TILT,
                         TILT_PRIOR,
                         PSI,
                         PSI_PRIOR,
                         SHIFTX,
                         SHIFTY,
                         SHIFTZ]

RELION_40_TOMO_LABELS = [TOMO_NAME,
                         SUBTOMO_NAME,
                         CTF_IMAGE,
                         COORD_X,
                         COORD_Y,
                         COORD_Z,
                         SHIFTX_ANGST,
                         SHIFTY_ANGST,
                         SHIFTZ_ANGST,
                         ROT,
                         TILT,
                         TILT_PRIOR,
                         PSI,
                         PSI_PRIOR,
                         CLASS_NUMBER]

# IMOD's eTomo files
NEWST_COM = 'newst.com'
TILT_COM = 'tilt.com'

# Data files
IN_TOMOS_STAR = 'inTomograms.star'
IN_TS_STAR = 'inTiltSeries.star'
OUT_TOMOS_STAR = 'tomograms.star'
IN_COORDS_STAR = 'inCoords.star'
IN_PARTICLES_STAR = 'inParticles.star'
OUT_PARTICLES_STAR = 'particles.star'
TRAJECTORIES_STAR = 'trajectories.star'
MANIFOLDS_STAR = 'manifolds.star'

REFINE_FSC_REF_STAR = '_model.star'
REFINE_STAR_FSC_TABLE = 'model_class_1'
REFINE_STAR_FSC_COLUMNS = ['rlnGoldStandardFsc']

FRAMES_DIR = 'frames'
MOTIONCORR_DIR = 'motioncorr'
TOMOGRAMS_DIR = 'tomograms'
POSTPROCESS_DIR = 'PostProcess'
FSC_REF_STAR = join(POSTPROCESS_DIR, 'postprocess.star')
POSTPROCESS_STAR_FSC_TABLE = 'fsc'
POSTPROCESS_STAR_FSC_COLUMNS = ['rlnFourierShellCorrelationCorrected',
                                'rlnFourierShellCorrelationUnmaskedMaps',
                                'rlnFourierShellCorrelationMaskedMaps',
                                'rlnCorrectedFourierShellCorrelationPhaseRandomizedMaskedMaps']

RE4_INDIV_GEN_FILES = [OUT_TOMOS_STAR, OUT_PARTICLES_STAR, TRAJECTORIES_STAR, MANIFOLDS_STAR, FSC_REF_STAR]
OPTIMISATION_SET_STAR = 'optimisation_set.star'
PSUBTOMOS_SQLITE = 'pseudosubtomograms%s.sqlite'


# Initial models
REC_PARTICLES_DIR = 'recParticles'
PSEUDO_SUBTOMOS_DIR = 'pseudoSubtomos'

# Refine - angular sampling
ANGULAR_SAMPLING_LIST = ['30', '15', '7.5', '3.7', '1.8', '0.9', '0.5',
                         '0.2', '0.1', '0.06', '0.03', '0.01', '0.007', '0.004']

# Symmetry description
jsonData = '[{"group": "Asymmetric", "notation": "C1", "origin": "User-defined", "orientation": "User-defined"},' \
           '{"group": "Cyclic", "notation": "C<n>", "origin": "On symm axis, Z user-defined", "orientation": "Symm axis on Z"},' \
           '{"group": "Dihedral", "notation": "D<n>", "origin": "Intersection of symm axes", "orientation": "principle symm axis on Z, 2-fold on X"},' \
           '{"group": "Tetrahedral", "notation": "T", "origin": "Intersection of symm axes", "orientation": "3-fold axis on Z (deviating from Heymann et al!)"},' \
           '{"group": "Octahedral", "notation": "O", "origin": "Intersection of symm axes", "orientation": "4-fold axes on X, Y, Z"},' \
           '{"group": "Icosahedral", "notation": "I<n>", "origin": "Intersection of symm axes", "orientation": "**"}]'
SYMMETRY_TABLE = json.loads(jsonData)

SYMMETRY_HELP_MSG = 'Symmetry libraries have been copied from XMIPP. As such, with the exception of tetrahedral ' \
                    'symmetry, they comply with ' \
                    'https://relion.readthedocs.io/en/latest/Reference/Bibliography.html#id23. Possible values ' \
                    '[notation label] are described below:\n\n' \
                    '%s' % json.dumps(SYMMETRY_TABLE, indent=1)


# Per particle per tilt box size values
BOX_SIZE_VALS = list(range(32, 512, 16))

# Aberration orders to be estimated in CTF refine
ABERRATION_ORDERS_EVEN = [4, 6, 8]
ABERRATION_ORDERS_ODD = [3, 5, 7]

# Star file comparer main messages
STAR_FILES_EQUAL = 'The introduced star files are equal.'
STAR_DIFF_SIZE = 'Different number of rows found:'
STAR_DIFF_LABELS = 'Different labels found:'
STAR_DIFF_VALUES = 'Different values found:'

# Tomograms star tables
GLOBAL_TABLE = 'global'

# Particles star tables
OPTICS_TABLE = 'optics'
PARTICLES_TABLE = 'particles'

# Post process files
POST_PROCESS_MRC = 'postprocess.mrc'
POST_PROCESS_MASKED_MRC = 'postprocess_masked.mrc'

# RELION5 ##############################################################################################################
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

alignedTsStarFields4TomoRec = [
    RLN_MICROGRAPH_MOVIENAME,
    RLN_TOMO_TILT_MOVIE_FRAME_COUNT,
    RLN_TOMO_NOMINAL_STAGE_TILT_ANGLE,
    RLN_TOMO_NOMINAL_TILT_AXIS_ANGLE,
    RLN_MICROGRAPH_PRE_EXPOSURE,
    # RLN_TOMO_NOMINAL_DEFOCUS,
    # RLN_CTF_POWER_SPECTRUM,
    RLN_MICROGRAPH_NAME_EVEN,
    RLN_MICROGRAPH_NAME_ODD,
    RLN_MICROGRAPH_NAME,
    RLN_MICROGRAPH_METADATA,
    RLN_ACCUM_MOTION_TOTAL,
    RLN_ACCUM_MOTION_EARLY,
    RLN_ACCUM_MOTION_LATE,
    # RLN_CTF_IMAGE,
    # RLN_DEFOCUS_U,
    # RLN_DEFOCUS_V,
    # RLN_CTF_ASTIGMATISM,
    # RLN_DEFOCUS_ANGLE,
    # RLN_CTF_FIGURE_OF_MERIT,
    # RLN_CTF_MAX_RESOLUTION,
    # RLN_CTF_ICE_RING_DENSITY,
    RLN_TOMO_X_TILT,
    RLN_TOMO_Y_TILT,
    RLN_TOMO_Z_ROT,
    RLN_TOMO_X_SHIFT_ANGST,
    RLN_TOMO_Y_SHIFT_ANGST,
    # RLN_CTF_SCALEFACTOR
]


# TiltImage extended fields for Relion
RELION_MIC_MOVIE_NAME = '_rlnMicrographMovieName'
RELION_MOT_CORR_IMG_NAME = '_rlnMicrographName'
RELION_NOMINAL_DEFOCUS = '_rlnTomoNominalDefocus'
RELION_MIC_METADATA = '_rlnMicrographMetadata'
RELION_ACCUM_DOSE_TOTAL = '_rlnAccumMotionTotal'
RELION_ACCUM_DOSE_EARLY = '_rlnAccumMotionEarly'
RELION_ACCUM_DOSE_LATER = '_rlnAccumMotionLate'

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

alignedTsSetStarFields4TomoRec = [
    RLN_TOMONAME,
    RLN_VOLTAGE,
    RLN_SPHERICALABERRATION,
    RLN_AMPLITUDECONTRAST,
    RLN_MICROGRAPHORIGINALPIXELSIZE,
    RLN_TOMOHAND,
    RLN_OPTICSGROUPNAME,
    RLN_TOMOTILT_SERIES_PIXEL_SIZE,
    RLN_TOMOTILT_SERIES_STAR_FILE,
    RLN_ETOMO_DIRECTIVE_FILE
]

# TILT-SERIES MOVIES METADATA ##########################################################################################
# Set file
tsMStarFields = [
    RLN_TOMONAME,
    RLN_TOMOTILT_SERIES_STAR_FILE,
    RLN_VOLTAGE,
    RLN_SPHERICALABERRATION,
    RLN_AMPLITUDECONTRAST,
    RLN_MICROGRAPHORIGINALPIXELSIZE,
    RLN_TOMOHAND,
    RLN_OPTICSGROUPNAME
]
# Each TSM file
tsMTsStarFields = [
    RLN_MICROGRAPH_MOVIENAME,
    RLN_TOMO_TILT_MOVIE_FRAME_COUNT,
    RLN_TOMO_NOMINAL_STAGE_TILT_ANGLE,
    RLN_TOMO_NOMINAL_TILT_AXIS_ANGLE,
    RLN_MICROGRAPH_PRE_EXPOSURE,
    RLN_TOMO_NOMINAL_DEFOCUS,
]

# COORDINATES METADATA #################################################################################################
GENERAL_TABLE = 'general'

RLN_ARE2DSTACKS = 'rlnTomoSubTomosAre2DStacks'

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

# Relion 5 coordinates attributes
R5_ROT_ATTRIB = '_rlnAngleRot'
R5_TILT_ATTRIB = '_rlnAngleTilt'
R5_PSI_ATTRIB = '_rlnAnglePsi'
R5_TILT_PRIOR_ATTRIB = '_rlnAngleTiltPrior'
R5_PSI_PRIO_ATTRIB = '_rlnAnglePsiPrior'
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

POSTPROCESS_STAR_FIELD = '_postprocessStar'
