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

from emtable import Table
from reliontomo import Plugin
from reliontomo.constants import TOMO_NAME_30, PARTICLES_TABLE, RELION_30_TOMO_LABELS, RELION_40_TOMO_LABELS, \
    GENERAL_TABLE, TOMO_PARTICLE_ID
from reliontomo.convert import convert40_tomo, convert30_tomo, convert50_tomo


def createWriterTomo(isPyseg=False, **kwargs):
    if isPyseg:
        writer = createWriterTomo30(starHeaders=RELION_30_TOMO_LABELS, **kwargs)
    elif Plugin.IS_GT50():
        writer = createWriterTomo50()
    else:
        writer = createWriterTomo40(starHeaders=RELION_40_TOMO_LABELS, **kwargs)
    return writer


def createWriterTomo30(**kwargs):
    """ Create a new Writer instance for Relion 3."""
    Writer = convert30_tomo.Writer
    return Writer(**kwargs)


def createWriterTomo40(**kwargs):
    """ Create a new Writer instance for Relion 4."""
    Writer = convert40_tomo.Writer
    return Writer(**kwargs)


def createWriterTomo50(**kwargs):
    """ Create a new Writer instance for Relion 5."""
    Writer = convert50_tomo.Writer
    return Writer(**kwargs)


def writeSetOfCoordinates(coordSet, starFile, whitelist, **kwargs):
    return createWriterTomo40().coordinates2Star(coordSet, starFile, whitelist, **kwargs)


def writeSetOfSubtomograms(particlesSet, starFile, **kwargs):
    """ Convenience function to write a SetOfSubtomograms as Relion metadata using a Writer."""
    writer = createWriterTomo(**kwargs)
    return writer.subtomograms2Star(particlesSet, starFile)


def writeSetOfPseudoSubtomograms(particlesSet, starFile, withPriors=False, **kwargs):
    return createWriterTomo40(**kwargs).pseudoSubtomograms2Star(particlesSet, starFile, withPriors=withPriors)


def writeSetOfTomograms(imgSet, starFile, whiteList=None, **kwargs):
    """ Convenience function to write a SetOfTomograms as Relion metadata using a Writer."""
    # There's no tomogram star file in Relion3
    return createWriterTomo40(**kwargs).tiltSeries2Star(imgSet, starFile, whiteList=whiteList, **kwargs)


def createReaderTomo(starFile, isRelion5, tableName=PARTICLES_TABLE, isCoordsStar=False,  **kwargs):
    dataTable = Table()
    try:
        dataTable.read(starFile, tableName=tableName)
    except:
        dataTable.read(starFile, tableName=None)

    labels = dataTable.getColumnNames()
    if TOMO_NAME_30 in labels:
        reader = convert30_tomo.Reader(starFile, dataTable, **kwargs)
        readerVersion = 3
    else:
        if isCoordsStar:
            if TOMO_PARTICLE_ID in labels:
                reader = convert40_tomo.Reader(starFile, dataTable, **kwargs)
                readerVersion = 4
            else:
                reader = convert50_tomo.Reader(starFile, dataTable, **kwargs)
                readerVersion = 5
        else:
            if isRelion5:
                reader = convert50_tomo.Reader(starFile, dataTable, **kwargs)
                readerVersion = 5
            else:
                reader = convert40_tomo.Reader(starFile, dataTable, **kwargs)
                readerVersion = 4
    return reader, readerVersion


def readSetOfPseudoSubtomograms(outputSet, isRelion5=True):
    """ Convenience function to write a SetOfPseudoSubtomograms as Relion metadata using a Reader."""
    # Subtomograms are represented in Relion 4 as Pseudosubtomograms
    reader, _ = createReaderTomo(outputSet.getParticlesStar(), isRelion5=isRelion5)
    return reader.starFile2PseudoSubtomograms(outputSet)


def readTsStarFile(inTs, outTs, starFile, outStackName, extraPath, isEvenOdd=False):
    dataTable = Table()
    dataTable.read(starFile, tableName=outTs.getTsId())
    reader = convert50_tomo.Reader(starFile, dataTable)
    return reader.starFile2Ts(inTs, outTs, outStackName, extraPath, isEvenOdd=isEvenOdd)
