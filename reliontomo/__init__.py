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

import pwem
import relion
from pwem import Config
from pyworkflow.utils import Environ
from reliontomo.constants import V3_0, RELIONTOMO_HOME, RELIONTOMO_DEFAULT, RELION, RELIONTOMO_CUDA_LIB

_logo = "relion_logo.png"
_references = ['Scheres2012a', 'Scheres2012b', 'Kimanius2016', 'Zivanov2018']
__version__ = '3.0.3'


class Plugin(relion.Plugin):
    _supportedVersions = [V3_0]
    _homeVar = RELIONTOMO_HOME
    _pathVars = [RELIONTOMO_HOME]

    @classmethod
    def _defineVariables(cls):
        cls._defineEmVar(RELIONTOMO_HOME, RELIONTOMO_DEFAULT)
        cls._defineVar(RELIONTOMO_CUDA_LIB, pwem.Config.CUDA_LIB)

    @classmethod
    def IS_30(cls):
        return RELIONTOMO_DEFAULT in cls.getVar(RELIONTOMO_HOME)

    @classmethod
    def getEnviron(cls):
        """ Setup the environment variables needed to launch Relion. """
        environ = Environ(os.environ)
        binPath = os.pathsep.join([cls.getHome('bin'), Config.MPI_BINDIR])
        libPath = os.pathsep.join([cls.getHome('lib'), cls.getHome('lib64'), pwem.Config.MPI_LIBDIR])

        if binPath not in environ['PATH']:
            environ.update({'PATH': binPath,
                            'LD_LIBRARY_PATH': libPath
                            }, position=Environ.BEGIN)

        # Get Relion CUDA library path if defined
        cudaLib = cls.getVar(RELIONTOMO_CUDA_LIB, pwem.Config.CUDA_LIB)
        environ.addLibrary(cudaLib)

        if 'RELION_MPI_LIB' in os.environ:
            environ.addLibrary(os.environ['RELION_MPI_LIB'])

        if 'RELION_MPI_BIN' in os.environ:
            environ.set('PATH', os.environ['RELION_MPI_BIN'],
                        position=Environ.BEGIN)
        return environ

    @classmethod
    def defineBinaries(cls, env):
        relion_commands = [('cmake -DGUI=OFF -DCMAKE_INSTALL_PREFIX=./ .', []),
                           ('make -j %d' % env.getProcessors(),
                            ['bin/relion_refine'])]

        env.addPackage(RELION, version=V3_0,
                       url='https://github.com/3dem/relion/archive/3.0.tar.gz',
                       commands=relion_commands,
                       updateCuda=True,
                       default=True)

