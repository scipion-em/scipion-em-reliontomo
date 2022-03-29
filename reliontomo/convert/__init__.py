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
from reliontomo.constants import TOMO_NAME_30, COORD_X, COORD_Y, COORD_Z, SUBTOMO_NAME, CTF_MISSING_WEDGE, ROT, TILT, \
    PSI, SHIFTX, SHIFTY, SHIFTZ, PIXEL_SIZE, MAGNIFICATION, TILT_PRIOR, SHIFTX_ANGST, PSI_PRIOR, SHIFTY_ANGST, \
    SHIFTZ_ANGST, CTF_IMAGE, TOMO_NAME, CLASS_NUMBER
from reliontomo.convert import convert40_tomo, convert30_tomo

PYSEG_SUBTOMO_LABELS = [TOMO_NAME_30,
                        COORD_X,
                        COORD_Y,
                        COORD_Z,
                        SUBTOMO_NAME,
                        CTF_MISSING_WEDGE,
                        ROT,
                        TILT,
                        PSI,
                        SHIFTX,
                        SHIFTY,
                        SHIFTZ]

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


def createWriterTomo(**kwargs):
    if Plugin.isRe40():
        if kwargs.get('isPyseg', None):
            Writer = createWriterTomo30(starHeaders=PYSEG_SUBTOMO_LABELS, **kwargs)
        else:
            Writer = createWriterTomo40(starHeaders=RELION_40_TOMO_LABELS, **kwargs)
    else:
        Writer = createWriterTomo30(starHeaders=RELION_30_TOMO_LABELS, **kwargs)
    return Writer


def createWriterTomo30(**kwargs):
    """ Create a new Writer instance for Relion 3."""
    Writer = convert30_tomo.Writer
    return Writer(**kwargs)


def createWriterTomo40(**kwargs):
    """ Create a new Writer instance for Relion 4."""
    Writer = convert40_tomo.Writer
    return Writer(**kwargs)


def writeSetOfSubtomograms(particlesSet, starFile, **kwargs):
    """ Convenience function to write a SetOfSubtomograms as Relion metadata using a Writer."""
    writer = createWriterTomo(**kwargs)
    return writer.subtomograms2Star(particlesSet, starFile)


def writeSetOfPseudoSubtomograms(particlesSet, starFile, **kwargs):
    return createWriterTomo40(**kwargs).pseudoSubtomograms2Star(particlesSet, starFile)


def writeSetOfTomograms(imgSet, starFile, **kwargs):
    """ Convenience function to write a SetOfTomograms as Relion metadata using a Writer."""
    # There's no tomogram star file in Relion3
    return createWriterTomo40(**kwargs).tiltSeries2Star(imgSet, starFile, **kwargs)


def createReaderTomo(starFile=None, **kwargs):
    if starFile:
        dataTable = Table()
        dataTable.read(starFile)
        labels = dataTable.getColumnNames()
        if TOMO_NAME_30 in labels:
            reader = convert30_tomo.Reader(**kwargs)
            isReader40 = False
        else:
            reader = convert40_tomo.Reader(**kwargs)
            isReader40 = True
        reader.read(starFile)
    else:
        if Plugin.isRe40():
            reader = convert40_tomo.Reader(**kwargs)
            isReader40 = True
        else:
            reader = convert30_tomo.Reader(**kwargs)
            isReader40 = False

    return reader, isReader40


def readSetOfPseudoSubtomograms(starFile, outputSet):
    """ Convenience function to write a SetOfPseudoSubtomograms as Relion metadata using a Reader."""
    # Subtomogras are represented in Relion 4 as Pseudosubtomograms
    reader, _ = createReaderTomo()
    return reader.starFile2PseudoSubtomograms(starFile, outputSet)
