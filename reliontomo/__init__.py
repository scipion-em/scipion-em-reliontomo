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
import relion
import pyworkflow.utils as pwutils
from reliontomo.contants import V3_0, RELIONTOMO_HOME

_logo = "relion_logo.png"
_references = ['Scheres2012a', 'Scheres2012b', 'Kimanius2016', 'Zivanov2018']
__version__ = '3.0.1'


class Plugin(relion.Plugin):
    _supportedVersions = [V3_0]
    _homeVar = RELIONTOMO_HOME

    @classmethod
    def defineBinaries(cls, env):
        pass

    @classmethod
    def IS_30(cls):
        environ = pwutils.Environ(os.environ)
        return 'relion-3.0' in environ[RELIONTOMO_HOME]
