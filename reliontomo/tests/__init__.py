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
OUTPUT_COORDS = 'outputCoordinates'
RE4_TOMO = 're4tomo'
EMD_10439 = 'emd_10439'


class DataSetRe4Tomo(Enum):
    eTomoDir = 'eTomo'
    alignments = 'eTomo/TS_43.xf'
    tiltSeries = 'eTomo/TS_43.mrc'
    mdoc = 'eTomo/TS_43.mrc.mdoc'
    tomogram = 'TS_43.mrc'
    coordinates = 'coords.star'


class DataSetEmd10439(Enum):
    tomoEmd10439 = 'tomograms/emd_10439.mrc',
    coords3dStarFile = 'importFromStarFiles/picking_001_parts.star',
    coords3dStarFileWithSRate = 'importFromStarFiles/picking_001_parts_with_sRate.star',
    subtomogramsStarFile = 'importFromStarFiles/class_ap_r_ali_k1_split.star',
    scipionSqlite3dCoords = 'importFromScipionSqlite/coordinates.sqlite'


DataSet(name=RE4_TOMO, folder=RE4_TOMO, files={el.name: el.value for el in DataSetRe4Tomo})
DataSet(name=EMD_10439, folder=EMD_10439, files={el.name: el.value[0] for el in DataSetEmd10439})
