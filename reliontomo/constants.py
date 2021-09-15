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

from reliontomo.utils import genSymmetryTable

RELION = 'relion'
V4_0 = '4.0'
RELIONTOMO_HOME = 'RELIONTOMO_HOME'
RELIONTOMO_CUDA_LIB = 'RELION_CUDA_LIB'
RELIONTOMO_DEFAULT_VERSION = V4_0
RELIONTOMO_DEFAULT = RELION + '-' + RELIONTOMO_DEFAULT_VERSION
V30_VALIDATION_MSG = 'This version of Reliontomo plugin requires Relion 3.0 binaries. ' \
                     'Please install it if necessary via scipion3 installb relion-3.0'


# Relion 4 tomogram star file fields
TOMO_NAME = 'rlnTomoName'
TILT_SERIES_NAME = 'rlnTomoTiltSeriesName'
CTFPLOTTER_FILE = 'rlnTomoImportCtfPlotterFile'
IMOD_DIR = 'rlnTomoImportImodDir'
FRACTIONAL_DOSE = 'rlnTomoImportFractionalDose'
ACQ_ORDER_FILE = 'rlnTomoImportOrderList'
CULLED_FILE = 'rlnTomoImportCulledFile'

# Relion 4 subtomogram star file fields
SUBTOMO_NAME = 'rlnImageName'
COORD_X = 'rlnCoordinateX'
COORD_Y = 'rlnCoordinateY'
COORD_Z = 'rlnCoordinateZ'
SHIFTX = 'rlnOriginX'
SHIFTY = 'rlnOriginY'
SHIFTZ = 'rlnOriginZ'
ROT = 'rlnAngleRot'
TILT = 'rlnAngleTilt'
PSI = 'rlnAnglePsi'
TILT_PRIOR = 'rlnAngleTiltPrior'
PSI_PRIOR = 'rlnAnglePsiPrior'
CLASS_NUMBER = 'rlnClassNumber'

# Relion pseudo-subtomograms
TOMO_PARTICLE_NAME = 'rlnTomoParticleName'
RANDOM_SUBSET = 'rlnRandomSubset'
OPTICS_GROUP = 'rlnOpticsGroup'
SHIFTX_ANGST = 'rlnOriginXAngst'
SHIFTY_ANGST = 'rlnOriginYAngst'
SHIFTZ_ANGST = 'rlnOriginZAngst'
CTF_IMAGE = 'rlnCtfImage'

FILE_NOT_FOUND = 'File not found'

# IMOD's eTomo files
NEWST_COM = 'newst.com'
TILT_COM = 'tilt.com'

# Data preparation
IN_TOMOS_STAR = 'inTomos.star'
OUT_TOMOS_STAR = 'outTomos.star'
IN_SUBTOMOS_STAR = 'inSubtomos.star'
OUT_SUBTOMOS_STAR = 'particles.star'

# Initial models
REC_PARTICLES_DIR = 'recParticles'
PSEUDO_SUBTOMOS_DIR = 'pseudoSubtomos'
INITIAL_MODEL = 'initial_model.mrc'

# Refine - angular sampling
ANGULAR_SAMPLING_LIST = ['30', '15', '7.5', '3.7', '1.8', '0.9', '0.5',
                         '0.2', '0.1', '0.06', '0.03', '0.01', '0.007', '0.004']

# Symmetry description
SYMMETRY_HELP_MSG = 'Symmetry libraries have been copied from XMIPP. As such, with the exception of tetrahedral ' \
                    'symmetry, they comply with ' \
                    'https://relion.readthedocs.io/en/latest/Reference/Bibliography.html#id23. Possible values ' \
                    '[notation label] are described below:\n\n' \
                    '%s' % json.dumps(genSymmetryTable(), indent=1)



