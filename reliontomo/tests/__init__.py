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

from pyworkflow.tests import DataSet

DataSet(name='reliontomo', folder='reliontomo',
        files={
               'refVol': 'caja72_job30_run_it025_class005.mrc',
               'mask': 'mascara_cilindro_caja72.mrc',
               'tiltseries': '*.mrcs',
               'tomograms': '*ali_bin1.mrc',
               'tlts': '*.tlt',
               'doseFiles': '*.txt',
               'dynamoTables': '*.tbl',
               'coordsFromStarDir': 'importFromStarFiles',
               'tomoEmd10439': 'importFromStarFiles/tomograms/emd_10439.mrc',
               'coords3dStarFile': 'importFromStarFiles/picking_001_parts.star',
        })
