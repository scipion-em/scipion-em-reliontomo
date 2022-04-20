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
import logging
import csv
from emtable import Table
from pwem import ALIGN_NONE, ALIGN_2D, ALIGN_PROJ
from pwem.convert.headers import fixVolume
from pwem.objects import Transform
from pyworkflow.object import Integer
from pyworkflow.utils import getParentFolder, yellowStr
from relion.convert import OpticsGroups
from reliontomo.constants import TOMO_NAME, TILT_SERIES_NAME, CTFPLOTTER_FILE, IMOD_DIR, FRACTIONAL_DOSE, \
    ACQ_ORDER_FILE, CULLED_FILE, SUBTOMO_NAME, COORD_X, COORD_Y, COORD_Z, ROT, TILT, PSI, CLASS_NUMBER, \
    TOMO_PARTICLE_NAME, RANDOM_SUBSET, OPTICS_GROUP, SHIFTX_ANGST, SHIFTY_ANGST, SHIFTZ_ANGST, CTF_IMAGE, \
    TOMO_PARTICLE_ID, MANIFOLD_INDEX, MRC, FILE_NOT_FOUND, PARTICLES_TABLE
import pwem.convert.transformations as tfs
import numpy as np
from os.path import join
from reliontomo.convert.convertBase import checkSubtomogramFormat, getTransformInfoFromCoordOrSubtomo, WriterTomo, \
    ReaderTomo, getTransformMatrixFromRow
from reliontomo.objects import PSubtomogram
from tomo.constants import BOTTOM_LEFT_CORNER
from tomo.objects import Coordinate3D

logger = logging.getLogger(__name__)


RELION_3D_COORD_ORIGIN = BOTTOM_LEFT_CORNER


class Writer(WriterTomo):
    """ Helper class to convert from Scipion SetOfTomograms and SetOfSubTomograms to star files ."""
    
    def __init__(self,  **kwargs):
        super().__init__(**kwargs)

    def tiltSeries2Star(self, tsSet, outStarFileName, prot=None, ctfPlotterParentDir=None, eTomoParentDir=None):
        tsTable = Table(columns=self._getTomogramStarFileLabels())
        for ts in tsSet:
            tsId = ts.getTsId()
            tsTable.addRow(tsId,                                                # _rlnTomoName #1
                           ts.getFirstItem().getFileName() + ':mrc',            # _rlnTomoTiltSeriesName #2
                           self._getCtfPlotterFile(tsId, ctfPlotterParentDir),  # _rlnTomoImportCtfPlotterFile #3
                           join(eTomoParentDir, tsId),                          # _rlnTomoImportImodDir #4
                           ts.getAcquisition().getDosePerFrame(),               # _rlnTomoImportFractionalDose #5
                           self._genOrderListFile(prot, ts),                    # _rlnTomoImportOrderList #6
                           self._genCulledFileName(prot, tsId)                  # _rlnTomoImportCulledFile #7
                           )

        # Write the STAR file
        tsTable.write(outStarFileName)

    def coordinates2Star(self, coordSet, subtomosStar, sRate=1, coordsScale=1):
        """Input coordsScale is used to scale the coordinates so they are expressed in bin 1, as expected by Relion 4"""
        tomoTable = Table(columns=self._getCoordinatesStarFileLabels())
        i = 0
        for coord in coordSet.iterCoordinates():
            angles, shifts = getTransformInfoFromCoordOrSubtomo(coord)
            # Add row to the table which will be used to generate the STAR file
            tomoTable.addRow(
                coord.getTomoId(),                                           # 1 _rlnTomoName
                coord.getObjId(),                                            # 2 _rlnTomoParticleId
                coord.getGroupId() if coord.getGroupId() else 1,             # 3 _rlnTomoManifoldIndex
                # coord in pix at scale of bin1
                coord.getX(RELION_3D_COORD_ORIGIN) * coordsScale,            # 4 _rlnCoordinateX
                coord.getY(RELION_3D_COORD_ORIGIN) * coordsScale,            # 5 _rlnCoordinateY
                coord.getZ(RELION_3D_COORD_ORIGIN) * coordsScale,            # 6 _rlnCoordinateZ
                # pix * Å/pix = [shifts in Å]
                shifts[0] * sRate,                                           # 7 _rlnOriginXAngst
                shifts[1] * sRate,                                           # 8 _rlnOriginYAngst
                shifts[2] * sRate,                                           # 9 _rlnOriginZAngst
                # Angles in degrees
                angles[0],                                                   # 10 _rlnAngleRot
                angles[1],                                                   # 11 _rlnAngleTilt
                angles[2],                                                   # 12 _rlnAnglePsi
                # Extended fields
                int(getattr(coord, '_classNumber', -1)),                     # 13_rlnClassNumber
                # Alternated 1 and 2 values
                int(getattr(coord, '_randomSubset', (i % 2) + 1)),           # 14 _rlnRandomSubset
            )
            i += 1
        # Write the STAR file
        tomoTable.write(subtomosStar)

    def pseudoSubtomograms2Star(self, pSubtomoSet, subtomosStar):
        sRate = pSubtomoSet.getSamplingRate()
        tomoTable = Table(columns=self._getPseudoSubtomogramStarFileLabels())

        # Write the STAR file
        optGroup = OpticsGroups.fromString(pSubtomoSet.getAcquisition().opticsGroupInfo.get())
        with open(subtomosStar, 'w') as f:
            optGroup.toStar(f)
            # Write header first
            partsWriter = Table.Writer(f)
            partsWriter.writeTableName('particles')
            partsWriter.writeHeader(tomoTable.getColumns())
            for subtomo in pSubtomoSet:
                coord = subtomo._coordinate
                angles, shifts = getTransformInfoFromCoordOrSubtomo(subtomo)
                # Add row to the table which will be used to generate the STAR file
                partsWriter.writeRowValues([
                    coord.getTomoId(),                                           # _rlnTomoName #1
                    subtomo.getObjId(),                                          # rlnTomoParticleId #2
                    coord.getGroupId(),                                          # rlnTomoManifoldIndex #3
                    # Coords in pixels
                    coord.getX(RELION_3D_COORD_ORIGIN),                          # _rlnCoordinateX #4
                    coord.getY(RELION_3D_COORD_ORIGIN),                          # _rlnCoordinateY #5
                    coord.getZ(RELION_3D_COORD_ORIGIN),                          # _rlnCoordinateZ #6
                    # pix * Å/pix = [shifts in Å]
                    shifts[0] * sRate,                                           # _rlnOriginXAngst #7
                    shifts[1] * sRate,                                           # _rlnOriginYAngst #8
                    shifts[2] * sRate,                                           # _rlnOriginZAngst #9
                    # Angles in degrees
                    angles[0],                                                   # _rlnAngleRot #10
                    angles[1],                                                   # _rlnAngleTilt #11
                    angles[2],                                                   # _rlnAnglePsi #12

                    subtomo.getClassId(),                                        # _rlnClassNumber #13
                    subtomo._randomSubset.get(),                                 # _rlnRandomSubset #14
                    subtomo._particleName.get(),                                 # _rlnTomoParticleName #15
                    subtomo._opticsGroupId.get(),                                # _rlnOpticsGroup #16
                    subtomo.getFileName().replace(':' + MRC, ''),                # _rlnImageName #17
                    subtomo._ctfImage.get()                                      # _rlnCtfImage #18
                ])
                
    def subtomograms2Star(self, subtomoSet, subtomosStar):
        logger.info("Writing relion4 star file (%s) from subtomograms."% subtomosStar)
        tomoTable = Table(columns=self.starHeaders)
        sRate = subtomoSet.getSamplingRate()
        extraPath = join(getParentFolder(subtomosStar), 'extra')
        for subtomo in subtomoSet.iterSubtomos():
            checkSubtomogramFormat(subtomo, extraPath)
            angles, shifts = getTransformInfoFromCoordOrSubtomo(subtomo)
            ctfFile = getattr(subtomo, '_ctfImage', None)
            if ctfFile:
                ctfFile = ctfFile.get()
            classNumber = subtomo.getClassId()
            rlnClassNumber = classNumber if classNumber else 1
            rlnTomoName = subtomo.getVolName()
            rlnImageName = subtomo.getFileName().replace(':' + MRC, '')
            rlnCtfImage = ctfFile if ctfFile else FILE_NOT_FOUND
            # Coords in pixels
            rlnCoordinateX = 0
            rlnCoordinateY = 0
            rlnCoordinateZ = 0

            if subtomo.hasCoordinate3D():
                rlnCoordinateX = subtomo.getCoordinate3D().getX(BOTTOM_LEFT_CORNER)
                rlnCoordinateY = subtomo.getCoordinate3D().getY(BOTTOM_LEFT_CORNER)
                rlnCoordinateZ = subtomo.getCoordinate3D().getZ(BOTTOM_LEFT_CORNER)

            rlnAngleRot = angles[0]
            rlnAngleTilt = angles[1]
            rlnAnglePsi = angles[2]
            # pix * Å/pix = [shifts in Å]
            rlnOriginX = shifts[0] * sRate
            rlnOriginY = shifts[1] * sRate
            rlnOriginZ = shifts[2] * sRate
            # Angles in degrees
            rlnTiltPrior = subtomo._tiltPriorAngle.get() if hasattr(subtomo, '_tiltPriorAngle') else rlnAngleTilt
            rlnPsiPrior = subtomo._psiPriorAngle.get() if hasattr(subtomo, '_psiPriorAngle') else rlnAnglePsi
            # Add row to the table which will be used to generate the STAR file
            fieldsToAdd = [rlnTomoName,
                           rlnImageName,
                           rlnCtfImage,
                           rlnCoordinateX,
                           rlnCoordinateY,
                           rlnCoordinateZ,
                           rlnOriginX,
                           rlnOriginY,
                           rlnOriginZ,
                           rlnAngleRot,
                           rlnAngleTilt,
                           rlnTiltPrior,
                           rlnAnglePsi,
                           rlnPsiPrior,
                           rlnClassNumber]

            tomoTable.addRow(*fieldsToAdd)

        # Write the STAR file
        tomoTable.write(subtomosStar)

    @staticmethod
    def _getTomogramStarFileLabels():
        return [
            TOMO_NAME,
            TILT_SERIES_NAME,
            CTFPLOTTER_FILE,
            IMOD_DIR,
            FRACTIONAL_DOSE,
            ACQ_ORDER_FILE,
            CULLED_FILE
        ]

    @staticmethod
    def _getCoordinatesStarFileLabels():
        return [
            TOMO_NAME,
            TOMO_PARTICLE_ID,
            MANIFOLD_INDEX,
            COORD_X,
            COORD_Y,
            COORD_Z,
            SHIFTX_ANGST,
            SHIFTY_ANGST,
            SHIFTZ_ANGST,
            ROT,
            TILT,
            PSI,
            CLASS_NUMBER,
            RANDOM_SUBSET
        ]

    @staticmethod
    def _getPseudoSubtomogramStarFileLabels():
        pSubtomosLabels = Writer._getCoordinatesStarFileLabels()
        pSubtomosLabels.extend([
            CLASS_NUMBER,
            RANDOM_SUBSET,
            TOMO_PARTICLE_NAME,
            OPTICS_GROUP,
            SUBTOMO_NAME,
            CTF_IMAGE
        ])
        return pSubtomosLabels

    @staticmethod
    def _getCtfPlotterFile(tsId, ctfPlotterDir):
        return join(ctfPlotterDir, tsId, tsId + '.defocus')

    @staticmethod
    def _genOrderListFile(prot, ts):
        """The order file expected by Relion is A 2-column, comma-separated file with the frame-order list
        of the tilt series, where the first column is the frame (image) number (starting at 1) and the second
        column is the tilt angle (in degrees).
        :param prot: current protocol object
        :param ts: TiltSeries object"""
        outputFilename = prot._getExtraPath(ts.getTsId() + '_order_list.csv')
        tiList = [ti.clone() for ti in ts]
        ind = np.argsort([ti.getAcquisitionOrder() for ti in tiList])  # Indices to get the data sorted by acqOrder
        with open(outputFilename, mode='w') as acqOrderFile:
            acqOrderFileWriter = csv.writer(acqOrderFile, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
            acqOrderList = [ti.getAcquisitionOrder() for ti in tiList]
            if min(acqOrderList) == 0:
                [acqOrderFileWriter.writerow([tiList[i].getAcquisitionOrder() + 1, tiList[i].getTiltAngle()]) for i in ind]
            else:
                [acqOrderFileWriter.writerow([tiList[i].getAcquisitionOrder(), tiList[i].getTiltAngle()]) for i in ind]

        return outputFilename

    @staticmethod
    def _genCulledFileName(prot, tsId):
        return prot._getExtraPath(tsId + '_culled.mrc:mrc')


class Reader(ReaderTomo):

    ALIGNMENT_LABELS = [
        SHIFTX_ANGST,
        SHIFTY_ANGST,
        SHIFTZ_ANGST,
        ROT,
        TILT,
        PSI,
    ]

    def __init__(self, starFile, dataTable, **kwargs):
        super().__init__(starFile, dataTable)
        self._shifts = np.zeros(3)
        self._angles = np.zeros(3)
        self._alignType = kwargs.get('alignType', ALIGN_NONE)
        self._pixelSize = kwargs.get('pixelSize', 1.0)

    @staticmethod
    def gen3dCoordFromStarRow(row, sRate, precedentIdDict, factor=1):
        coordinate3d = None
        tomoId = row.get(TOMO_NAME)
        vol = precedentIdDict.get(tomoId, None)
        if vol:
            coordinate3d = Coordinate3D()
            x = row.get(COORD_X, 0)
            y = row.get(COORD_Y, 0)
            z = row.get(COORD_Z, 0)
            coordinate3d.setVolume(vol)
            coordinate3d.setX(float(x)*factor, RELION_3D_COORD_ORIGIN)
            coordinate3d.setY(float(y)*factor, RELION_3D_COORD_ORIGIN)
            coordinate3d.setZ(float(z)*factor, RELION_3D_COORD_ORIGIN)
            coordinate3d.setMatrix(getTransformMatrixFromRow(row, sRate=sRate))
            coordinate3d.setGroupId(row.get(MANIFOLD_INDEX, 1))
            # Extended fields
            coordinate3d._classNumber = Integer(row.get(CLASS_NUMBER, -1))
            coordinate3d._randomSubset = Integer(row.get(RANDOM_SUBSET, 1))

        return coordinate3d, tomoId

    def starFile2Coords3D(self, coordsSet, precedentsSet, scaleFactor):
        precedentIdDict = {}
        for tomo in precedentsSet:
            precedentIdDict[tomo.getTsId()] = tomo.clone()

        nonMatchingTomoIds = ''
        for row in self.dataTable:
            # Consider that there can be coordinates in the star file that does not belong to any of the tomograms
            # introduced
            coord, tomoId = self.gen3dCoordFromStarRow(row, coordsSet.getSamplingRate(),
                                                       precedentIdDict, factor=scaleFactor)
            if coord:
                coordsSet.append(coord)
            else:
                if tomoId not in nonMatchingTomoIds:
                    nonMatchingTomoIds += '%s ' % tomoId

        if nonMatchingTomoIds:
            print(yellowStr('The star file contains coordinates that belong to tomograms not present '
                            'in the introduced set of tomograms: %s' % nonMatchingTomoIds))

    def starFile2PseudoSubtomograms(self, outputSet):
        samplingRate = outputSet.getSamplingRate()
        listOfFilesToFixVolume = []

        for counter, row in enumerate(self.dataTable):
            particleFile = row.get(SUBTOMO_NAME)
            ctfFile = row.get(CTF_IMAGE)
            psubtomo = PSubtomogram(fileName=particleFile,
                                    samplingRate=samplingRate,
                                    ctfFile=ctfFile,
                                    tsId=row.get(TOMO_NAME),
                                    classId=row.get(CLASS_NUMBER, -1))
            # Add the files to the list of files whose header has to be corrected to be interpreted as volumes
            listOfFilesToFixVolume.append(particleFile)
            listOfFilesToFixVolume.append(ctfFile)
            # Add current pseudosubtomogram to the output set
            outputSet.append(psubtomo)

        # Fix volume headers
        fixVolume(listOfFilesToFixVolume)

    def setParticleTransform(self, particle, row, sRate):
        """ Set the transform values from the row. """

        if (self._alignType == ALIGN_NONE) or not row.hasAnyColumn(self.ALIGNMENT_LABELS):
            self.setParticleTransform = self.__setParticleTransformNone
        else:
            # Ensure the Transform object exists
            self._angles = np.zeros(3)
            self._shifts = np.zeros(3)

            particle.setTransform(Transform())

            if self._alignType == ALIGN_2D:
                self.setParticleTransform = self.__setParticleTransform2D
            elif self._alignType == ALIGN_PROJ:
                self.setParticleTransform = self.__setParticleTransformProj
            else:
                raise TypeError("Unexpected alignment type: %s"
                                % self._alignType)

        # Call again the modified function
        self.setParticleTransform(particle, row, sRate)

    @staticmethod
    def __setParticleTransformNone(particle, row, sRate):
        particle.setTransform(None)

    def __setParticleTransform2D(self, particle, row, sRate):
        angles = self._angles
        shifts = self._shifts

        shifts[0] = float(row.get(SHIFTX_ANGST / sRate, 0))
        shifts[1] = float(row.get(SHIFTY_ANGST / sRate, 0))
        angles[2] = float(row.get(PSI, 0))
        radAngles = -np.deg2rad(angles)
        M = tfs.euler_matrix(radAngles[0], radAngles[1], radAngles[2], 'szyz')
        M[:3, 3] = shifts[:3]
        particle.getTransform().setMatrix(M)

    def __setParticleTransformProj(self, particle, row, sRate):
        angles = self._angles
        shifts = self._shifts

        shifts[0] = float(row.get(SHIFTX_ANGST, 0)) / sRate
        shifts[1] = float(row.get(SHIFTY_ANGST, 0)) / sRate
        shifts[2] = float(row.get(SHIFTZ_ANGST, 0)) / sRate

        angles[0] = float(row.get(ROT, 0))
        angles[1] = float(row.get(TILT, 0))
        angles[2] = float(row.get(PSI, 0))

        radAngles = -np.deg2rad(angles)
        M = tfs.euler_matrix(radAngles[0], radAngles[1], radAngles[2], 'szyz')
        M[:3, 3] = -shifts[:3]
        M = np.linalg.inv(M)
        particle.getTransform().setMatrix(M)

