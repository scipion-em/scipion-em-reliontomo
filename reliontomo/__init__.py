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
import pwem
import relion
from reliontomo.constants import RELIONTOMO_HOME, RELIONTOMO_DEFAULT, RELION, RELIONTOMO_CUDA_LIB, V4_0

_logo = "relion_logo.png"
_references = ['Scheres2012a', 'Scheres2012b', 'Kimanius2016', 'Zivanov2018']
__version__ = '4.0.0'


class Plugin(relion.Plugin):
    _supportedVersions = [V4_0]
    _homeVar = RELIONTOMO_HOME
    _pathVars = [RELIONTOMO_HOME]
    _isReconstringParticleFromSubtomos = False

    @classmethod
    def _defineVariables(cls):
        cls._defineEmVar(RELIONTOMO_HOME, 'relion-%s' % V4_0)
        cls._defineVar(RELIONTOMO_CUDA_LIB, pwem.Config.CUDA_LIB)

    @staticmethod
    def isRe40():
        return True if Plugin.getHome().endswith(V4_0) else False

    @classmethod
    def runRelionTomo(cls, protocol, program, args, cwd=None, numberOfMpi=1):
        """ Run Relion command from a given protocol. """
        protocol.runJob(program, args, cwd=cwd, env=cls.getEnviron(), numberOfMpi=numberOfMpi)

    @classmethod
    def isReconstringParticleFromSubtomos(cls):
        return cls._isReconstringParticleFromSubtomos
