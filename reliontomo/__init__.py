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
_logo = "relion_logo.png"
_references = ['Scheres2012a', 'Scheres2012b', 'Kimanius2016', 'Zivanov2018']
__version__ = '3.2.1'

try:
    import pwem
    from reliontomo.constants import RELIONTOMO_HOME, RELIONTOMO_DEFAULT, RELION, RELIONTOMO_CUDA_LIB, V4_0
    import relion
    from relion import V3_1

    class Plugin(relion.Plugin):
        _supportedVersions = [V4_0]
        _homeVar = RELIONTOMO_HOME
        _pathVars = [RELIONTOMO_HOME]

        @classmethod
        def _defineVariables(cls):
            cls._defineEmVar(RELIONTOMO_HOME, 'relion-%s' % V4_0)
            cls._defineVar(RELIONTOMO_CUDA_LIB, pwem.Config.CUDA_LIB)

        @staticmethod
        def isRe40():
            """We only allow one version older than 4.0, which is 3.1"""
            return False if Plugin.getHome().endswith(V3_1) else True

        @classmethod
        def runRelionTomo(cls, protocol, program, args, cwd=None, numberOfMpi=1):
            """ Run Relion command from a given protocol. """
            protocol.runJob(program, args, cwd=cwd, env=cls.getEnviron(), numberOfMpi=numberOfMpi)

        @classmethod
        def defineBinaries(cls, env):
            pass

except Exception as e:
    pass
