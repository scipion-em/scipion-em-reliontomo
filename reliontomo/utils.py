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
import numpy as np
import emtable
from collections import OrderedDict
from os.path import isabs, join

def getProgram(program, nMpi):
    """ Get the program name depending on the MPI use or not."""
    if nMpi > 1:
        program += '_mpi'
    return program


def manageDims(fileName, z, n):
    if fileName.endswith('.mrc') or fileName.endswith('.map'):
        if z == 1 and n != 1:
            zDim = n
        else:
            zDim = z
    else:
        zDim = z

    return zDim


def getAbsPath(starFilePath, tomoFile):
    """If the paths of the files pointed from a star file are relative, they'll be referred to the path of the
    star file. This method is used to consider that case."""
    if isabs(tomoFile):
        return tomoFile
    else:
        return join(starFilePath, tomoFile)


def _gen2LevelBaseName(fullFileName):
    """This function generates a new basename composed of the original basename preceded by the parent folder
    followed by '_'. This will reduce the probabilities of non-uniqueness of the name when creating links."""
    parts = fullFileName.split('/')
    if len(parts) == 1:
        return fullFileName
    else:
        return parts[-2] + '_' + parts[-1]


def isPseudoSubtomogram(subtomo):
    return hasattr(subtomo, '_ctfImage')


def genEnumParamDict(keyList):
    return OrderedDict((key, val) for key, val in zip(keyList, range(len(keyList))))


def generateProjections(particlesFilePath, tomogramsFilePath):
    mdFileName = '%s@%s' % ('particles', particlesFilePath)
    table = emtable.Table(fileName=particlesFilePath)
    particles = []
    tomoNames = []

    for row in table.iterRows(mdFileName):
        if not row.get('rlnTomoName') in tomoNames:
            tomoNames.append(row.get('rlnTomoName'))
        particles.append([row.get('rlnTomoName'),
                          np.array([row.get('rlnCoordinateX'),
                                    row.get('rlnCoordinateY'),
                                    row.get('rlnCoordinateZ'), 1])])

    tomograms = {}
    for tomoName in tomoNames:
        mdFileName = '%s@%s' % (tomoName, tomogramsFilePath)
        table = emtable.Table(fileName=tomogramsFilePath)
        projections = []
        for row in table.iterRows(mdFileName):
            projections.append(np.array([eval(row.get('rlnTomoProjX')),
                                         eval(row.get('rlnTomoProjY')),
                                         eval(row.get('rlnTomoProjZ')),
                                         eval(row.get('rlnTomoProjW'))]))
        tomograms[tomoName] = projections

    return projectParticles(particles, tomograms)


def projectParticles(particles, tomograms):
    projections = []

    for partId, particle in enumerate(particles):
        tomoName = particle[0]
        tomoProjections = tomograms[tomoName]
        for tiltId, tomoProjection in enumerate(tomoProjections):
            multproj = tomoProjection.dot(particle[1])
            projections.append([tomoName,  tiltId, partId, multproj[0], multproj[1]])

    return projections
