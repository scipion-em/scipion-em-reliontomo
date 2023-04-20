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
IN_TOMOS_STAR = 'inTomos.star'
OUT_TOMOS_STAR = 'tomograms.star'
IN_COORDS_STAR = 'inCoords.star'
IN_PARTICLES_STAR = 'inSubtomos.star'
OUT_PARTICLES_STAR = 'particles.star'
TRAJECTORIES_STAR = 'trajectories.star'
MANIFOLDS_STAR = 'manifolds.star'
POSTPROCESS_DIR = 'PostProcess'
FSC_REF_STAR = join(POSTPROCESS_DIR, 'postprocess.star')
RE4_INDIV_GEN_FILES = [OUT_TOMOS_STAR, OUT_PARTICLES_STAR, TRAJECTORIES_STAR, MANIFOLDS_STAR, FSC_REF_STAR]
OPTIMISATION_SET_STAR = 'optimisation_set.star'
PSUBTOMOS_SQLITE = 'pseudosubtomograms%s.sqlite'


# Initial models
REC_PARTICLES_DIR = 'recParticles'
PSEUDO_SUBTOMOS_DIR = 'pseudoSubtomos'
INITIAL_MODEL = 'initial_model.mrc'

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
