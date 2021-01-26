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
V30_VALIDATION_MSG = 'This version of Reliontomo plugin requires Relion 3.0. Please install it if necessary and set ' \
                     'variable RELIONTOMO_HOME in scipion.conf and pass it when launching scipion3:' \
                     ' ./scipion3 --config [path_to_scipion.conf].'
RELIONTOMO_HOME = 'RELIONTOMO_HOME'
RELIONTOMO_DEFAULT_VERSION = V3_0
RELIONTOMO_DEFAULT = RELION + '-' + RELIONTOMO_DEFAULT_VERSION

