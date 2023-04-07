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
logger = logging.getLogger(__name__)
import csv

from emtable import Table
from pwem import ALIGN_NONE
from pwem.convert.headers import fixVolume
from pwem.objects import Transform
from pyworkflow.object import Integer, Float
from pyworkflow.utils import getParentFolder, yellowStr, createLink, makePath
from relion.convert import OpticsGroups
from reliontomo.constants import *
import numpy as np
from os.path import join, basename
from reliontomo.convert.convertBase import (checkSubtomogramFormat,
                                            getTransformInfoFromCoordOrSubtomo,
                                            WriterTomo, ReaderTomo, getTransformMatrixFromRow)
from reliontomo.objects import RelionPSubtomogram
from tomo.constants import BOTTOM_LEFT_CORNER, TR_SCIPION, TR_RELION, SCIPION
from tomo.objects import Coordinate3D, SubTomogram, TomoAcquisition

logger = logging.getLogger(__name__)


class Writer(WriterTomo):
    """ Helper class to convert from Scipion SetOfTomograms and SetOfSubTomograms to star files ."""

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    def tiltSeries2Star(self, tsSet, outStarFileName, prot=None, ctfPlotterParentDir=None, eTomoParentDir=None,
                        whiteList=None):
        """Writes the needed tomograms.star file for relion prepare
        :param tsSet: set of tilt series
        :param outStarFileName: name of the star file that will be generated
        :param prot: protocol object
        :param ctfPlotterParentDir: ctf plotter results parent directory path
        :param eTomoParentDir: eTomo results parent directory path
        :param whiteList list with tilseriesIds to be filter by"""

        tsTable = Table(columns=self._getTomogramStarFileLabels())
        for ts in tsSet:
            tsId = ts.getTsId()

            if whiteList is None or tsId in whiteList:
                tsTable.addRow(tsId,  # _rlnTomoName #1
                               ts.getFirstItem().getFileName(),  # _rlnTomoTiltSeriesName #2
                               self._getCtfPlotterFile(tsId, ctfPlotterParentDir),  # _rlnTomoImportCtfPlotterFile #3
                               join(eTomoParentDir, tsId),  # _rlnTomoImportImodDir #4
                               ts.getAcquisition().getDosePerFrame(),  # _rlnTomoImportFractionalDose #5
                               self._genOrderListFile(prot, ts),  # _rlnTomoImportOrderList #6
                               self._genCulledFileName(prot, tsId)  # _rlnTomoImportCulledFile #7
                               )
            else:
                logger.info("Tilt series %s excluded." % tsId)
        # Write the STAR file
        tsTable.write(outStarFileName)

    def coordinates2Star(self, coordSet, subtomosStar, sRate=1, coordsScale=1):
        """Input coordsScale is used to scale the coordinates so they are expressed in bin 1, as expected by Relion 4"""
        tomoTable = Table(columns=self._getCoordinatesStarFileLabels())
        i = 0
        for coord in coordSet.iterCoordinates():
            angles, shifts = getTransformInfoFromCoordOrSubtomo(coord, coordSet.getSamplingRate())
            # Add row to the table which will be used to generate the STAR file
            tomoTable.addRow(
                coord.getTomoId(),  # 1 _rlnTomoName
                coord.getObjId(),  # 2 _rlnTomoParticleId
                coord.getGroupId() if coord.getGroupId() else 1,  # 3 _rlnTomoManifoldIndex
                # coord in pix at scale of bin1
                coord.getX(RELION_3D_COORD_ORIGIN) * coordsScale,  # 4 _rlnCoordinateX
                coord.getY(RELION_3D_COORD_ORIGIN) * coordsScale,  # 5 _rlnCoordinateY
                coord.getZ(RELION_3D_COORD_ORIGIN) * coordsScale,  # 6 _rlnCoordinateZ
                # pix * Å/pix = [shifts in Å]
                shifts[0], #* sRate,  # 7 _rlnOriginXAngst
                shifts[1], #* sRate,  # 8 _rlnOriginYAngst
                shifts[2], #* sRate,  # 9 _rlnOriginZAngst
                # Angles in degrees
                angles[0],  # 10 _rlnAngleRot
                angles[1],  # 11 _rlnAngleTilt
                angles[2],  # 12 _rlnAnglePsi
                # Extended fields
                int(getattr(coord, '_classNumber', -1)),  # 13_rlnClassNumber
                # Alternated 1 and 2 values
                int(getattr(coord, '_randomSubset', (i % 2) + 1)),  # 14 _rlnRandomSubset
                int(coord.getX(SCIPION)),  # 15 _sciXCoord
                int(coord.getY(SCIPION)),  # 16 _sciYCoord
                int(coord.getZ(SCIPION))   # 17 _sciZCoord
            )
            i += 1
        # Write the STAR file
        tomoTable.write(subtomosStar)

    def pseudoSubtomograms2Star(self, pSubtomoSet, outStar, withPriors=False):

        logger.info("Generating particles file (%s) from pseudosubtomogram set." % outStar)
        sRate = pSubtomoSet.getSamplingRate()
        hasCoords = pSubtomoSet.getFirstItem().hasCoordinate3D()
        tomoTable = Table(columns=self._getPseudoSubtomogramStarFileLabels(hasCoords, withPriors=withPriors))

        # Write the STAR file
        optGroup = OpticsGroups.fromString(pSubtomoSet.getAcquisition().opticsGroupInfo.get())
        with open(outStar, 'w') as f:
            optGroup.toStar(f)
            # Write header first
            partsWriter = Table.Writer(f)
            partsWriter.writeTableName(PARTICLES_TABLE)
            partsWriter.writeHeader(tomoTable.getColumns())
            for pSubtomo in pSubtomoSet:
                angles, shifts = getTransformInfoFromCoordOrSubtomo(pSubtomo, pSubtomo.getSamplingRate(), convention=TR_RELION)
                pSubtomoFile = pSubtomo.getFileName()
                pSubtomoFile = pSubtomoFile.replace(':' + MRC, '') if pSubtomoFile else FILE_NOT_FOUND
                pSubtomoCtfFile = pSubtomo.getCtfFile() if pSubtomo.getCtfFile() else FILE_NOT_FOUND

                # Add row to the table which will be used to generate the STAR file
                rowsValues = [
                    pSubtomo.getTsId(),  # _rlnTomoName #1
                    pSubtomo.getObjId(),  # rlnTomoParticleId #2
                    pSubtomo.getManifoldIndex(),  # rlnTomoManifoldIndex #3
                    # Coords in pixels
                    pSubtomo.getX(),  # _rlnCoordinateX #4
                    pSubtomo.getY(),  # _rlnCoordinateY #5
                    pSubtomo.getZ(),  # _rlnCoordinateZ #6
                    # pix * Å/pix = [shifts in Å]
                    shifts[0], #* sRate,  # _rlnOriginXAngst #7
                    shifts[1], #* sRate,  # _rlnOriginYAngst #8
                    shifts[2], #* sRate,  # _rlnOriginZAngst #9
                    # Angles in degrees
                    angles[0],  # _rlnAngleRot #10
                    angles[1],  # _rlnAngleTilt #11
                    angles[2],  # _rlnAnglePsi #12

                    pSubtomo.getClassId(),  # _rlnClassNumber #13
                    pSubtomo.getRdnSubset(),  # _rlnRandomSubset #14
                ]

                if withPriors:
                    rowsValues += [angles[0], angles[1], angles[2]]
                if hasCoords:
                    rowsValues += [pSubtomo.getCoordinate3D().getX(SCIPION),  # _sciXCoord #15
                                   pSubtomo.getCoordinate3D().getY(SCIPION),  # _sciYCoord #16
                                   pSubtomo.getCoordinate3D().getZ(SCIPION)]  # _sciZCoord #17

                rowsValues += [pSubtomo.getRe4ParticleName(),  # _rlnTomoParticleName #18
                               pSubtomo.getOpticsGroupId(),  # _rlnOpticsGroup #19
                               pSubtomoFile,  # _rlnImageName #20
                               pSubtomoCtfFile]  # _rlnCtfImage #21

                partsWriter.writeRowValues(rowsValues)

    def subtomograms2Star(self, subtomoSet, subtomosStar):
        logger.info("Writing relion4 star file (%s) from subtomograms." % subtomosStar)
        tomoTable = Table(columns=self.starHeaders)
        sRate = subtomoSet.getSamplingRate()
        extraPath = join(getParentFolder(subtomosStar), 'extra')
        for subtomo in subtomoSet.iterSubtomos():
            checkSubtomogramFormat(subtomo, extraPath)
            angles, shifts = getTransformInfoFromCoordOrSubtomo(subtomo, subtomo.getSamplingRate())
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
            rlnOriginX = shifts[0] #* sRate
            rlnOriginY = shifts[1] #* sRate
            rlnOriginZ = shifts[2] #* sRate
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
    def _getCoordinatesStarFileLabels(hasCoords=True, withPriors=False):
        starFileLabels = [
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

        if withPriors:
            starFileLabels+=[ROT_PRIOR, TILT_PRIOR, PSI_PRIOR]

        if hasCoords:
            starFileLabels += [SCIPION_COORD_X, SCIPION_COORD_Y, SCIPION_COORD_Z]

        return starFileLabels

    @staticmethod
    def _getPseudoSubtomogramStarFileLabels(hasCoords=True, withPriors=False):
        pSubtomosLabels = Writer._getCoordinatesStarFileLabels(hasCoords, withPriors=withPriors)
        pSubtomosLabels.extend([
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
                [acqOrderFileWriter.writerow([tiList[i].getAcquisitionOrder() + 1, tiList[i].getTiltAngle()]) for i in
                 ind]
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
            coordinate3d.setX(float(x) * factor, RELION_3D_COORD_ORIGIN)
            coordinate3d.setY(float(y) * factor, RELION_3D_COORD_ORIGIN)
            coordinate3d.setZ(float(z) * factor, RELION_3D_COORD_ORIGIN)
            coordinate3d.setMatrix(getTransformMatrixFromRow(row, sRate=sRate), convention=TR_RELION)
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
            logger.info(yellowStr('The star file contains coordinates that belong to tomograms not present '
                            'in the introduced set of tomograms: %s' % nonMatchingTomoIds))

    def starFile2PseudoSubtomograms(self, outputSet):
        sRate = outputSet.getSamplingRate()
        listOfFilesToFixVolume = []
        for counter, row in enumerate(self.dataTable):
            t = Transform()
            particleFile = row.get(SUBTOMO_NAME, None)
            ctfFile = row.get(CTF_IMAGE, None)
            psubtomo = RelionPSubtomogram(fileName=particleFile,
                                          samplingRate=sRate,
                                          ctfFile=ctfFile,
                                          tsId=row.get(TOMO_NAME),
                                          classId=row.get(CLASS_NUMBER, -1),
                                          x=row.get(COORD_X),
                                          y=row.get(COORD_Y),
                                          z=row.get(COORD_Z),
                                          rdnSubset=row.get(RANDOM_SUBSET, counter % 2),  # 1 and 2 alt. by default
                                          re4ParticleName=row.get(TOMO_PARTICLE_NAME),
                                          opticsGroupId=row.get(OPTICS_GROUP, 1),
                                          manifoldIndex=row.get(MANIFOLD_INDEX, 1 if counter % 2 else -1),
                                          )  # 1 and -1

            # Keeping particle id
            psubtomo.setObjId(row.get(TOMO_PARTICLE_ID))

            # Set the coordinate3D
            if row.get(SCIPION_COORD_X) is not None:  # Assume that the coordinates exists
                sciCoord = Coordinate3D()
                sciCoord.setX(row.get(SCIPION_COORD_X), SCIPION)
                sciCoord.setY(row.get(SCIPION_COORD_Y), SCIPION)
                sciCoord.setZ(row.get(SCIPION_COORD_Z), SCIPION)
                sciCoord.setTomoId(row.get(TOMO_NAME))
                psubtomo.setCoordinate3D(sciCoord)

            # Set the transformation matrix
            t.setMatrix(getTransformMatrixFromRow(row, sRate=sRate))
            psubtomo.setTransform(t)

            # Add the files to the list of files whose header has to be corrected to be interpreted as volumes
            if particleFile is not None and particleFile != FILE_NOT_FOUND:
                listOfFilesToFixVolume.append(particleFile)
            if ctfFile is not None and ctfFile!=FILE_NOT_FOUND:
                listOfFilesToFixVolume.append(ctfFile)
            # Add current pseudosubtomogram to the output set
            outputSet.append(psubtomo)

        # Keep the number of particles to compare sizes in case of subset
        outputSet.setNReParticles(self.dataTable.size())
        # Fix volume headers
        if listOfFilesToFixVolume:
            fixVolume(listOfFilesToFixVolume)

    def starFile2SubtomogramsImport(self, subtomoSet, coordSet, linkedSubtomosDir, starFilePath):
        samplingRate = subtomoSet.getSamplingRate()
        precedentsDict = {tomo.getTsId(): tomo for tomo in coordSet.getPrecedents()}
        for row, coordinate3d in zip(self.dataTable, coordSet):
            subtomo = SubTomogram()
            transform = Transform()
            origin = Transform()

            # Files
            tomoId = row.get(TOMO_NAME, FILE_NOT_FOUND)
            tomoName = precedentsDict[tomoId].getFileName()
            subtomoName = row.get(SUBTOMO_NAME, FILE_NOT_FOUND)
            tomoFolder = join(linkedSubtomosDir, tomoId)
            makePath(tomoFolder)
            linkedSubtomoName = join(tomoFolder, basename(subtomoName))
            createLink(subtomoName, linkedSubtomoName)  # Link the subtomos to the extra folder

            # Subtomograms
            tiltPrior = row.get(TILT_PRIOR, 0)
            psiPrior = row.get(PSI_PRIOR, 0)
            transform.setMatrix(coordinate3d.getMatrix(convention=TR_RELION))

            subtomo.setVolName(tomoName)
            subtomo.setFileName(linkedSubtomoName)
            subtomo.setCoordinate3D(coordinate3d)
            subtomo.setVolId(coordinate3d.getVolId())
            subtomo.setTransform(transform)
            subtomo.setAcquisition(TomoAcquisition())
            subtomo.setClassId(row.get(CLASS_NUMBER, 0))
            subtomo.setSamplingRate(samplingRate)
            subtomo._tiltPriorAngle = Float(tiltPrior)
            subtomo._psiPriorAngle = Float(psiPrior)
            subtomo.setOrigin(origin)

            # Add current subtomogram to the output set
            subtomoSet.append(subtomo)

        # Set the set of coordinates which corresponds to the current set of subtomograms
        subtomoSet.setCoordinates3D(coordSet)

