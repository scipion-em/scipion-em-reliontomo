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
import argparse

from reliontomo.constants import STAR_COMPARE_ENTRY_POINT
from reliontomo.objects import StarFileComparer


def main():
    parser = argparse.ArgumentParser(prog=STAR_COMPARE_ENTRY_POINT,
                                     usage="Reliontomo program to compare two star files")
    parser.add_argument('--star1', type=str, help='First star file to be compared.')
    parser.add_argument('--star2', type=str, help='Second star file to be compared.')
    parser.add_argument('--tableList', nargs='+', default=['particles'],
                        help='List of tables desired to be compared in both star files. Default is "particles".')
    parser.add_argument('--excludeLabelList', nargs='+', default=[],
                        help='list of labels to be excluded in the comparison.')

    args = parser.parse_args()
    star1 = args.star1
    star2 = args.star2
    tableList = args.tableList
    excludeLabelList = args.excludeLabelList

    for tableName in tableList:
        sfc = StarFileComparer(star1, star2, tableName)
        print(sfc.compare(excludeLabelsList=excludeLabelList))


if __name__ == '__main__':
    main()
