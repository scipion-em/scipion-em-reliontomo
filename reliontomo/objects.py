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
from enum import Enum
from os.path import exists, join

from emtable import Table

from pyworkflow.object import String, Integer, Float
from relion.convert import OpticsGroups
from reliontomo.constants import (OPT_TOMOS_STAR, OPT_PARTICLES_STAR,
                                  OPT_TRAJECTORIES_STAR, OPT_MANIFOLDS_STAR,
                                  OPT_FSC_STAR, OUT_TOMOS_STAR, OUT_PARTICLES_STAR,
                                  TRAJECTORIES_STAR, MANIFOLDS_STAR, FSC_REF_STAR,
                                  STAR_DIFF_SIZE, STAR_DIFF_LABELS, STAR_DIFF_VALUES,
                                  STAR_FILES_EQUAL, PSUBTOMOS_SQLITE, OPTICS_TABLE)
from tomo.constants import SCIPION, TR_SCIPION
from tomo.objects import SetOfSubTomograms, SubTomogram


class EnumRe4GenFilesProps(Enum):
    _tomograms = OUT_TOMOS_STAR
    _particles = OUT_PARTICLES_STAR
    _trajectories = TRAJECTORIES_STAR
    _manifolds = MANIFOLDS_STAR
    _referenceFsc = FSC_REF_STAR


class RelionPSubtomogram(SubTomogram):

    TS_ID_ATTRIBUTE = '_volName'

    def __init__(self, fileName=None, samplingRate=None, ctfFile=None, tsId=None, classId=None,
                 x=None, y=None, z=None, rdnSubset=None, re4ParticleName=None,
                 opticsGroupId=1, manifoldIndex=None, logLikeliCont=None, maxValProbDist=None,
                 noSignifSamples=None, **kwargs):
        super().__init__(**kwargs)
        self.setFileName(fileName)
        self.setSamplingRate(samplingRate)
        self.setClassId(classId)
        self._volName.set(tsId)
        self._ctfFile = String(ctfFile)
        self._x = Float(x)
        self._y = Float(y)
        self._z = Float(z)
        self._rdnSubset = Integer(rdnSubset)
        self._re4ParticleName = String(re4ParticleName)
        self._opticsGroupId = Integer(opticsGroupId)
        self._manifoldIndex = Integer(manifoldIndex)
        self._logLikeliContribution = Float(logLikeliCont)
        self._maxValueProbDistribution = Float(maxValProbDist)
        self._nrOfSignificantSamples = Integer(noSignifSamples)

    def getCtfFile(self):
        return self._ctfFile.get()

    def getTsId(self):
        return self.getVolName()

    def getBoxSize(self):
        return self._boxSize.get()

    def getX(self):
        return self._x.get()

    def getY(self):
        return self._y.get()

    def getZ(self):
        return self._z.get()

    def getCoords(self):
        return self.getX(), self.getY(), self.getZ()

    def getRdnSubset(self):
        return self._rdnSubset.get()

    def getRe4ParticleName(self):
        return self._re4ParticleName.get()

    def getOpticsGroupId(self):
        return self._opticsGroupId.get()

    def getManifoldIndex(self):
        return self._manifoldIndex.get()

    def setCtfFile(self, val):
        self._ctfFile.set(val)

    def setTsId(self, val):
        self._volName.set(val)

    def setTransform(self, newTransform, convention=TR_SCIPION):
        super().setTransform(newTransform, convention=convention)

    def setBoxSize(self, val):
        self._boxSize = val

    def setX(self, val):
        self._x = val

    def setY(self, val):
        self._y = val

    def setZ(self, val):
        self._z = val

    def setCoords(self, x, y, z):
        self.setX(x)
        self.setY(y)
        self.setZ(z)

    def setRdnSubset(self, val):
        self._rdnSubset.set(val)

    def setReParticleName(self, val):
        self._re4ParticleName.set(val)

    def setOpticsGroupId(self, val):
        self._opticsGroupId.set(val)

    def setManifoldIndex(self, val):
        self._manifoldIndex.set(val)

    # def copyInfo(self, other):
    #     self.copyAttributes(other, '_samplingRate', '_tsId', '_rdnSubset',
    #                         '_re4ParticleName', '_opticsGroupId', '_boxSize')


class RelionSetOfPseudoSubtomograms(SetOfSubTomograms):
    """ Set to persist relion's metadata files and particles.
    Approach: heep always the generated optimization_set.star file. Additionally,
    Keep the particles.star file is apply (e.g.: apply operation)

    Run relion with optimization star file and additionally, particles if is not none (empty it when apply)

    For subsets... we could generate always a new particles star file, either based on this set
    (when we fixed the shifts conversion  problem) or alternatively read the relion starfile and generate
    a new one filtered out base on  this set ids? (warning, ids will fail if there are joins!!)
    """
    ITEM_TYPE = RelionPSubtomogram

    def __init__(self, optimSetStar=None, relionBinning=None, tsSamplingRate=None, boxSize=24,
                 nReParticles=0, **kwargs):
        super().__init__(**kwargs)
        self._filesMaster = String() # Optimisation set file path
        self._boxSize = Integer(boxSize)
        self._tomograms = String() # Path to tomograms star file, usually the one prepared by relion at the beginning
        self._particles = String() # Path to particles.star file where particles metadata
        self._trajectories = String() # Path to trajectories file
        self._manifolds = String() # Path to manifolds file
        self._referenceFsc = String() # FSC file
        self._relionBinning = Float(relionBinning) # Binning of this set
        self._tsSamplingRate = Float(tsSamplingRate) # Sampling rate of the tilt series
        self._nReParticles = Integer(nReParticles) # number of relion particles in the particles star file

        if optimSetStar:
            self.filesMaster = optimSetStar

    # def __str__(self):
    #     strRep = ''
    #     if self.getNumParticles():
    #         strRep += '%i items, ' % self.getNumParticles()
    #     if self.getRelionBinning():
    #         strRep += 'binning %.1f' % self.getRelionBinning()
    #     return self.getClassName() + ' %s' % ('(%s)' % strRep if strRep else '')

    @property
    def filesMaster(self):
        return self._filesMaster.get()

    @filesMaster.setter
    def filesMaster(self, optimSetStar):
        try:
            self._readOptimSetStar(optimSetStar)
            self._filesMaster.set(optimSetStar)
        except FileNotFoundError:
            raise FileNotFoundError('Unable to find file %s' % optimSetStar)
        except TypeError:
            raise TypeError('No optimisation set star file was provided.')

    def _readOptimSetStar(self, optimSetStar):
        dataTable = Table()
        dataTable.read(optimSetStar)
        rowObj = dataTable[0]  # This file only contains the different filenames related to the current STA step
        self._tomograms.set(rowObj.get(OPT_TOMOS_STAR, None))
        self._particles.set(rowObj.get(OPT_PARTICLES_STAR, None))
        self._trajectories.set(rowObj.get(OPT_TRAJECTORIES_STAR, None))
        self._manifolds.set(rowObj.get(OPT_MANIFOLDS_STAR, None))
        self._referenceFsc.set(rowObj.get(OPT_FSC_STAR, None))
        # Read optimisation set and fill the corresponding attribute
        og = OpticsGroups(Table(fileName=self._particles.get(), tableName=OPTICS_TABLE))
        acq = self.getAcquisition()
        acq.opticsGroupInfo = String(og.toString())
        self.setAcquisition(acq)

    def updateGenFiles(self, extraPath):
        """Some protocols don't generate the optimisation_set.star file. In that case, the input Object which
        represents it will be copied and, after that, this method will be used to update the corresponding
        attributes with the generated files."""
        for p in EnumRe4GenFilesProps:
            currentFile = join(extraPath, p.value)
            if exists(currentFile):
                setattr(self, p.name, String(currentFile))

    def getTomogramsStar(self):
        return self._tomograms.get()

    def getParticles(self):
        return self._particles.get()

    def getTrajectories(self):
        return self._trajectories.get()

    def getManifolds(self):
        return self._manifolds.get()

    def getReferenceFsc(self):
        return self._referenceFsc.get()

    def getRelionBinning(self):
        return self._relionBinning.get()

    def getTsSamplingRate(self):
        return self._tsSamplingRate.get()

    def getBoxSize(self):
        return self._boxSize.get()

    def getCurrentSamplingRate(self):
        return self.getTsSamplingRate() * self.getRelionBinning()

    def getNReParticles(self):
        return self._nReParticles.get()

    def setNReParticles(self, val):
        self._nReParticles.set(val)

    def setRelionBinning(self, val):
        self._relionBinning.set(val)

    def setTsSamplingRate(self, val):
        self._tsSamplingRate.set(val)

    def setBoxSize(self, val):
        self._boxSize.set(val)

    def copyInfo(self, other):
        self.copyAttributes(other,'_filesMaster', '_tomograms', '_particles', '_trajectories', '_manifolds', '_referenceFsc',
                            '_relionBinning', '_tsSamplingRate', '_samplingRate', '_boxSize', '_nReParticles',
                            '_coordsPointer')
        self._acquisition.copyInfo(other._acquisition)
        # self._relionMd = relionMd if relionMd else relionTomoMetadata
    def _samplingRateStr(self):
        return "%0.2f Å/px" % self.getCurrentSamplingRate()
        
    # def iterPSubtomos(self, ts=None, orderBy='id'):
    #     if ts is None:
    #         tsId = None
    #     elif isinstance(ts, str):
    #         tsId = ts
    #     elif isinstance(ts, TiltSeries):
    #         tsId = ts.getTsId()
    #     else:
    #         raise Exception('Invalid input tilt series of type %s' % type(ts))
    #
    #     # Iterate over all psubtomos if tomoId is None,
    #     # otherwise use tomoId to filter the where selection
    #     psubtomoWhere = '1' if tsId is None else '_tomoId=%s' % tsId
    #     for psubtomo in self.iterItems(where=psubtomoWhere, orderBy=orderBy):
    #         yield psubtomo
    #
    # def getRelionMd(self):
    #     return self._relionMd
    #
    # def setRelionMd(self, reMd):
    #     self._relionMd = reMd
    #
    # def copyInfo(self, other):
    #     self.copyAttributes(other, '_samplingRate', '_relionMd')
    #     self._acquisition.copyInfo(other._acquisition)


def createSetOfRelionPSubtomograms(protocolPath, optimSetStar, coordsPointer, template=PSUBTOMOS_SQLITE, tsSamplingRate=1,
                                   relionBinning=1, boxSize=24, nReParticles=0):
    """ Creates the RelionSetOfSubtomograms from the input arguments

    :param protocolPath: Path of the protocol where to create the sqlite
    :param optimSetStar: optimization set star file. This file is a small collection of files(images and metadata to pass to future jobs
    :param coordsPointer: Pointer to the set of coordinates
    :param tsSamplingRate: Sampling rate of the tilt series
    :param relionBinning: binning of the set
    :param boxSize: Box size of the set
    :param nReParticles: Number of particles in relion's particles star file

    """
    psubtomoSet = RelionSetOfPseudoSubtomograms.create(protocolPath, template=template)
    psubtomoSet.filesMaster = optimSetStar
    psubtomoSet.setTsSamplingRate(tsSamplingRate)
    psubtomoSet.setRelionBinning(relionBinning)
    psubtomoSet.setSamplingRate(psubtomoSet.getCurrentSamplingRate())
    psubtomoSet.setBoxSize(boxSize)
    psubtomoSet.setNReParticles(nReParticles)
    psubtomoSet.setCoordinates3D(coordsPointer)

    # Clone acquisition from tomograms
    acquisition = coordsPointer.get().getPrecedents().getAcquisition()
    newAcquisition = acquisition.clone()
    psubtomoSet.setAcquisition(newAcquisition)
    return psubtomoSet



class StarFileComparer:

    def __init__(self, starFile1, starFile2, table2Read):
        self._table2Read = table2Read
        self.dataTable1 = starFile1
        self.dataTable2 = starFile2
        self._labels = None

    def compare(self, excludeLabelsList=None):
        msg = ''
        msg += self.compareSize()
        if not msg:
            msg += self.compareLabels(excludeLabelsList=excludeLabelsList)
        if not msg:
            msg += self.compareValues(excludeLabelsList=excludeLabelsList)
        tableNameMsg = 'Table [%s]\n' % self._table2Read
        return tableNameMsg + msg if msg else tableNameMsg + STAR_FILES_EQUAL

    def compareSize(self):
        msg = ''
        mRows1 = len(self.dataTable1)
        nRows2 = len(self.dataTable2)
        if mRows1 != nRows2:
            msg = '\n- %s %i != %i' % (STAR_DIFF_SIZE, mRows1, nRows2)
        return msg

    def compareLabels(self, excludeLabelsList=None):
        msg = ''
        labels1 = self.dataTable1.getColumnNames()
        labels2 = self.dataTable2.getColumnNames()
        if excludeLabelsList:
            for label in excludeLabelsList:
                if label in labels1:
                    labels1.remove(label)
                if label in labels2:
                    labels2.remove(label)

        if labels1 == labels2:
            self._labels = labels1
        else:
            msg = '\n- %s' \
                  '\n\tLABELS 1: %s' \
                  '\n\tLABELS 2: %s' \
                  '\n\tDIFF LABELS: %s' % \
                  (STAR_DIFF_LABELS, list2str(labels1), list2str(labels2), list2str(set(labels1) ^ set(labels2)))
        return msg

    def compareValues(self, excludeLabelsList=None, tolerance=0.1):
        msg = ''
        labels = self._labels if not excludeLabelsList else self._updateLabelsList(self._labels, excludeLabelsList)
        counter = 1
        for row1, row2 in zip(self.dataTable1, self.dataTable2):
            rowMsg = ''
            for label in labels:
                val1 = row1.get(label)
                val2 = row2.get(label)
                try:
                    # Numeric case --> apply tolerance
                    val1 = float(val1)
                    val2 = float(val2)
                    if abs(float(val1) - float(val2)) > tolerance:
                        rowMsg += '\n\t\tLABEL = %s, %s != %s' % (label, val1, val2)
                except ValueError:
                    # Not numeric case  --> check equality
                    if val1 != val2:
                        rowMsg += '\n\t\tLABEL = %s, %s != %s' % (label, val1, val2)
            if rowMsg:
                msg += '\n\tROW %i%s' % (counter, rowMsg)
            counter += 1

        if msg:
            msg = STAR_DIFF_VALUES + msg
        return msg

    @property
    def dataTable1(self):
        return self._dataTable1

    @dataTable1.setter
    def dataTable1(self, value):
        self._dataTable1 = self._readStarFile(value)

    @property
    def dataTable2(self):
        return self._dataTable2

    @dataTable2.setter
    def dataTable2(self, value):
        self._dataTable2 = self._readStarFile(value)

    def _readStarFile(self, starFile):
        try:
            dataTable = Table()
            dataTable.read(starFile, tableName=self._table2Read)
            return dataTable

        except FileNotFoundError:
            raise FileNotFoundError('Unable to find file %s' % starFile)

        except TypeError:
            raise TypeError('A string containing he file and path of a star file was expected.')

    @staticmethod
    def _updateLabelsList(labelsList, excludeLabelsList):
        for label in excludeLabelsList:
            if label in labelsList:
                labelsList.remove(label)
        return labelsList


def list2str(inList):
    return ' '.join([str(label) for label in inList])



