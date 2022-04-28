# **************************************************************************
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
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************
from enum import Enum
from pyworkflow.tests import DataSet

OUTPUT_TOMOS = 'outputTomograms'
RE4_TOMO = 're4tomo'
EMD_10439 = 'emd_10439'


class DataSetRe4Tomo(Enum):
    eTomoDir = 'eTomo'
    alignments = '*/TS_*.xf'
    tiltSeries = '*/TS_*.mrc'
    mdocs = '*/*.mdoc'
    coordinates = 'coords.star'
    starFilesComparerDir = 'starFilesComparer'
    generalTestingDir = 'starFilesComparer/generalTesting'
    okStarFile = 'starFilesComparer/generalTesting/particles.star'
    diffSizeStar = 'starFilesComparer/generalTesting/particlesDiffSize5.star'
    diffLabelsStar = 'starFilesComparer/generalTesting/particlesDiffLabels.star'
    diffValuesStar = 'starFilesComparer/generalTesting/particlesDiffValues.star'
    relionStarDir = 'starFilesComparer/re4NativeStars'
    prepareTomosStarRelion = 'starFilesComparer/re4NativeStars/tomograms.star'
    prepareParticlesStarRelion = 'starFilesComparer/re4NativeStars/particles.star'
    makePSubtomosStarRelion = 'starFilesComparer/re4NativeStars/psubtomos.star'
    zShiftedRelion = 'starFilesComparer/re4NativeStars/native_z2_75.star'
    scipionStarDir = 'starFilesComparer/re4ScipionStars'
    prepareTomosStarScipion = 'starFilesComparer/re4ScipionStars/tomograms.star'
    preparePartcilesStarScipion = 'starFilesComparer/re4ScipionStars/particles.star'
    makePSubtomosStarScipion = 'starFilesComparer/re4ScipionStars/psubtomos.star'
    zShiftedScipion = 'starFilesComparer/re4ScipionStars/scipion_z2_75.star'


class DataSetEmd10439(Enum):
    tomoEmd10439 = 'tomograms/emd_10439.mrc',
    coords3dStarFile = 'importFromStarFiles/picking_001_parts.star',
    coords3dStarFileWithSRate = 'importFromStarFiles/picking_001_parts_with_sRate.star',
    subtomogramsStarFile = 'importFromStarFiles/class_ap_r_ali_k1_split.star',
    scipionSqlite3dCoords = 'importFromScipionSqlite/coordinates.sqlite'


DataSet(name=RE4_TOMO, folder=RE4_TOMO, files={el.name: el.value for el in DataSetRe4Tomo})
DataSet(name=EMD_10439, folder=EMD_10439, files={el.name: el.value[0] for el in DataSetEmd10439})
