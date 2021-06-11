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

RELION = 'relion'
V3_0 = '3.0'
RELIONTOMO_HOME = 'RELIONTOMO_HOME'
RELIONTOMO_CUDA_LIB = 'RELION_CUDA_LIB'
RELIONTOMO_DEFAULT_VERSION = V3_0
RELIONTOMO_DEFAULT = RELION + '-' + RELIONTOMO_DEFAULT_VERSION
V30_VALIDATION_MSG = 'This version of Reliontomo plugin requires Relion 3.0 binaries. ' \
                     'Please install it if necessary via scipion3 installb %s.' % RELIONTOMO_DEFAULT

