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
# *  e-mail address 'scipion-users@lists.sourceforge.net'
# *
# **************************************************************************
import glob

from pyworkflow import BETA
from tomo.protocols import ProtTomoBase
from reliontomo import Plugin
from relion.convert import Table
from os.path import join, abspath
import pyworkflow.utils.path as pwutils
from pwem.protocols import EMProtocol
from pyworkflow.protocol import params, StringParam, EnumParam, String, IntParam
import numpy as np
from operator import sub
import math
from pwem.emlib.image import ImageHandler
from reliontomo.constants import V30_VALIDATION_MSG

CTFDIRBASENAME = 'Ctf3D'
CTF3D_PER_VOLUME = 0
CTF3D_PER_SUBVOLUME = 1
# Dictionary keys for set of tilt series common params
SRATE = 'sRate'
VOLTAGE = 'voltage'
SPHAB = 'sphAb'
AMPCN = 'ampCn'
# CTF files extensions
CTFSTAR = '.star'
CTFMRC = '.mrc'


class ProtRelionEstimateCTF3D(EMProtocol, ProtTomoBase):
    """ Generates the CTF star and MRC files needed by relion for the CTF3D.
    """
    _label = 'CTF 3D estimation'
    _devStatus = BETA
    _outputClassName = 'SetOfCoordinates3D'

    def __init__(self, **args):
        EMProtocol.__init__(self, **args)
        self.bFactor = 4  # Info from Relion wiki
        self.tsExpandedList = []
        self.initialized = False
        self.ctfMRCFileList = []
        self._doseFromMdoc = None

    # --------------------------- DEFINE param functions --------------------------------------------
    def _defineParams(self, form):
        form.addSection(label='Inputs')
        form.addParam('inputCoordinates', params.PointerParam,
                      label='3D coordinates',
                      important=True,
                      pointerClass='SetOfCoordinates3D',
                      help='Select a set of subtomogram coordinates.')
        form.addParam('inputSetCTFTomoSeries', params.PointerParam,
                      label='CTF tomo series',
                      important=True,
                      pointerClass='SetOfCTFTomoSeries',
                      help='Select a set of CTF tomo series.')
        form.addParam('doseFilesPath', params.PathParam,
                      label="Dose files directory\n(only if not importing from mdoc)",
                      allowsNull=True,
                      help="Not necessary if the tilt series or tilt series movies contains the dose data. It's "
                           "the case when the they're imported using the mdoc files. Root directory of the dose "
                           "files for the tilt series.")
        form.addParam('filesPattern', StringParam,
                      label='Pattern',
                      allowsNull=True,
                      help="Pattern of the files to be imported.\n\n"
                           "The pattern can contain standard wildcards such as\n"
                           "*, ?, etc, or special ones like ### to mark some\n"
                           "digits in the filename as ID.\n\n"
                           "NOTE: wildcards and special characters "
                           "('*', '?', '#', ':', '%') cannot appear in the "
                           "actual path.")
        form.addParam('boxSize', IntParam,
                      label='Box Size',
                      important=True,
                      allowsNull=False,
                      help='Perform a 3D reconstruction from 2D CTF-images, with the given size in pixels')
        group = form.addGroup('CTF 3D Estimation Mode')
        group.addParam('ctf3dMode', EnumParam,
                       choices=self._getImportChoices(),
                       default=CTF3D_PER_VOLUME,
                       label='Choose CTF 3D estimation type',
                       help='CTF 3D can be estimated per volume (faster, usable in first iterations '
                            'of the processing procedure) or per subvolume (slower, used for the refinement).')
        form.addParallelSection(threads=3, mpi=1)

    # --------------------------- INSERT steps functions --------------------------------------------
    def _insertAllSteps(self):
        self._initialize()
        # Insert the steps
        writeDeps = self._insertFunctionStep("writeStarCtf3DStep")
        recFeps = self._insertFunctionStep("reconstructCtf3DStep", prerequisites=[writeDeps])
        self._insertFunctionStep('createOutputStep', prerequisites=[recFeps])

    # --------------------------- STEPS functions --------------------------------------------
    def writeStarCtf3DStep(self):
        tsCounter = 0
        sRate = self.tsSet.getSamplingRate()
        voltage = self.tsSet.getAcquisition().getVoltage()
        sphAb = self.tsSet.getAcquisition().getSphericalAberration()
        ampCn = self.tsSet.getAcquisition().getAmplitudeContrast()
        setTsInfo = {SRATE: sRate,
                     VOLTAGE: voltage,
                     SPHAB: sphAb,
                     AMPCN: ampCn}
        tsList = self.tsExpandedList if self.tsExpandedList else self.tsSet
        for ts in tsList:
            if self.EstimationMode == CTF3D_PER_VOLUME:
                self._estimateCTF3DPerVolume(ts, setTsInfo, tsCounter)
            else:
                self._estimateCTF3DPerSubvolume(ts, setTsInfo, tsCounter)
            tsCounter += 1

    def reconstructCtf3DStep(self):
        sRate = self.tsSet.getSamplingRate()
        boxSize = self.boxSize.get()
        program = "relion_reconstruct" if self.numberOfMpi == 1 else "relion_reconstruct_mpi"

        for tsExt in self.tsExpandedList:
            coordCounter = 0
            for ctfStarFile, ctfMRCFile in zip(tsExt.getCTFStarList(), tsExt.getCTFMRCList()):
                param = {"sampling": sRate,
                         "ctfStar": abspath(ctfStarFile),
                         "ctf3D": abspath(ctfMRCFile),
                         "boxSize": boxSize
                         }

                args = " --i %(ctfStar)s --o %(ctf3D)s --reconstruct_ctf %(boxSize)d --angpix %(sampling)f"
                self.runJob(program, args % param) #, env=Plugin.getEnviron())
                coordCounter += 1

    def createOutputStep(self):
        out_coords = self._createSetOfCoordinates3D(self.coordSet)  # Create an empty set of micrographs
        # Copy all the info of the inputs, then the mrc ctf star file attribute will added
        out_coords.copyInfo(self.coordSet)
        coordCounter = 0
        if self.EstimationMode == CTF3D_PER_VOLUME:
            for tsExp in self.tsExpandedList:
                coords = tsExp.getCoords()
                ctfMrc = tsExp.getCTFMRCList()[0]  # Only one element was created for each TS in per volume case
                for coord in coords:
                    coord.setObjId(coordCounter + 1)
                    coord._3dcftMrcFile = String(ctfMrc)
                    out_coords.append(coord)
                    coordCounter += 1
        else:
            for coord, ctfMrc in zip(self.coordSet, self.ctfMRCFileList):
                coord._3dcftMrcFile = String(ctfMrc)
                out_coords.append(coord)

        self._defineOutputs(outputCoordinates=out_coords)
        self._defineTransformRelation(self.inputCoordinates, out_coords)

    # --------------------------- INFO functions --------------------------------------------
    def _validate(self):
        """ Should be overriden in subclasses to
        return summary message for NORMAL EXECUTION.
        """
        validateMsgs = self._initialize()
        if not Plugin.IS_30():
            validateMsgs.append(V30_VALIDATION_MSG)
        return validateMsgs

    def _summary(self):
        summary = []

        return summary

    def _methods(self):
        return []

    # --------------------------- UTILS functions ---------------------------------------------------
    def _initialize(self):
        validateMsgs = []
        if not self.initialized:
            self.coordSet = self.inputCoordinates.get()
            self.ctfTomoSet = self.inputSetCTFTomoSeries.get()
            self.tsSet = self.ctfTomoSet.getSetOfTiltSeries()
            self.EstimationMode = self.ctf3dMode.get()

            # Check if the dose data is currently known (TS were imported from mdoc)
            self._doseFromMdoc = self._hasDosePerFrame()
            if self._doseFromMdoc:
                self._genTSExp()
            else:
                # Assign each dose file to a tilt series
                doseFilesNoOk = self._getDoseFiles()
                if doseFilesNoOk:
                    validateMsgs.append(doseFilesNoOk)
                self.initialized = True
        return validateMsgs

    def _hasDosePerFrame(self):
        # It's assumed that if the first tilt image of the first tilt series has dose per frame, all the rest of the
        # tilt series of the set will have that data, too
        return True if self.tsSet.getFirstItem().getFirstItem().getAcquisition().getDosePerFrame() else False

    def _getDoseFiles(self):
        path = self.doseFilesPath.get('').strip()
        pattern = self.filesPattern.get('').strip()
        wholePattern = join(path, pattern) if pattern else path
        matches = glob.glob(wholePattern)
        if matches:
            nTs = len(self.tsSet)
            nMatches = len(matches)
            if nTs == nMatches:
                self._assignDoseFilesToTS(matches)
            else:
                return ("The number of dose files must be equal to the number of tilt series contained in the "
                        "introduced set of tilt series.\nnDoseFiles [{}] != nTS [{}]").format(nTs, nMatches)
        else:
            return "There are no files matching the pattern %s" % pattern

    def _assignDoseFilesToTS(self, matches):
        msg = ''
        nonMatchingTS = []
        remTS = len(self.tsSet)
        remDose = len(matches)
        tsList = [ts.clone(ignoreAttrs=[]) for ts in self.tsSet]
        ctfSeriesList = [ctfSeries.clone(ignoreAttrs=[]) for ctfSeries in self.ctfTomoSet]
        tomoList = [tomo.clone() for tomo in self.coordSet.getPrecedents()]

        tomoList.sort(key=self._sortTomoNames)
        tsList.sort(key=self._sortIds)
        ctfSeriesList.sort(key=self._sortIds)
        matches.sort()
        counter = 0
        for ctfs, ts in zip(ctfSeriesList, tsList):
            # This clone command was used to pass the value by value instead of by reference, because
            # all the elements of list self.tsExpandedList were overwritten on each iteration of this loop
            tsId = ts.getTsId().replace('TS_', '')
            remTS -= 1
            for doseFile in matches:
                doseBaseName = pwutils.removeBaseExt(doseFile).replace('_ExpDose', '')
                if tsId in doseBaseName or doseBaseName in tsId:
                    # Get the corresponding subtomograms coordinates
                    coordList = [coord.clone() for coord in self.coordSet.iterCoordinates(volume=tomoList[counter])]
                    # Add to the TS Expanded list
                    self.tsExpandedList.append(ExtendedTS(ts, ctfs, DoseFile(doseFile), coordList))
                    matches.remove(doseFile)
                    remDose -= 1
                    break

            if remTS != remDose:
                nonMatchingTS.append('\n' + tsId)

            counter += 1

        if not self.tsExpandedList:
            msg += '\nNo matching dose files were found'
        elif nonMatchingTS:
            msg += "No matching dose file was found for the following TS:%s" % ''.join(i for i in nonMatchingTS)

    def _genTSExp(self):
        tsList = [ts.clone(ignoreAttrs=[]) for ts in self.tsSet]
        ctfSeriesList = [ctfSeries.clone(ignoreAttrs=[]) for ctfSeries in self.ctfTomoSet]
        tomoList = [tomo.clone() for tomo in self.coordSet.getPrecedents()]

        tomoList.sort(key=self._sortTomoNames)
        tsList.sort(key=self._sortIds)
        ctfSeriesList.sort(key=self._sortIds)
        counter = 0
        for ctfs, ts in zip(ctfSeriesList, tsList):
            # This clone command was used to pass the value by value instead of by reference, because
            # all the elements of list self.tsExpandedList were overwritten on each iteration of this loop
            # Get the corresponding subtomograms coordinates
            coordList = [coord.clone() for coord in self.coordSet.iterCoordinates(volume=tomoList[counter])]
            # Add to the TS Expanded list
            self.tsExpandedList.append(ExtendedTS(ts, ctfs, None, coordList))
            counter += 1

    @staticmethod
    def _sortTomoNames(tomoList):
        return sorted([tomo for tomo in tomoList.getFileName()])

    @staticmethod
    def _sortIds(obj):
        return obj.getTsId()

    def _getProgram(self, program='relion_preprocess'):
        """ Get the program name depending on the MPI use or not. """
        if self.numberOfMpi > 1:
            program += '_mpi'
        return program

    def getTSPath(self, tsFn):
        tomoBaseDir = pwutils.removeBaseExt(tsFn)
        pwutils.makePath(self._getExtraPath(tomoBaseDir))
        return self._getExtraPath(tomoBaseDir)

    def _getCtfFile(self, tsId, coordCounter=None, fileExt=CTFSTAR, ctfMode=CTF3D_PER_VOLUME):
        ctfDir = join(self.getTSPath(tsId), CTFDIRBASENAME)
        pwutils.makePath(ctfDir)
        if ctfMode == CTF3D_PER_VOLUME:
            baseCtfFn = pwutils.removeBaseExt(tsId) + "_ctf"
        else:
            baseCtfFn = pwutils.removeBaseExt(tsId) + "_ctf_%06d" % coordCounter
        return join(ctfDir, baseCtfFn + fileExt)

    @staticmethod
    # Get the Dose for the corresponding tilt angle via the lowest difference between the tilt angles
    # contained in the dose file and the tilt angle contained as an attribute for the current tilt image
    def _getCurrentDose(tiltAngleDegs, tiltList, doseList):
        diffList = list(map(lambda x: abs(sub(x, tiltAngleDegs)), tiltList))
        indMin = diffList.index(min(diffList))
        return doseList[indMin]

    @ staticmethod
    def _createTable():
        # Headers for Relion's CTS star file
        return Table(columns=['rlnDefocusU',
                              'rlnVoltage',
                              'rlnSphericalAberration',
                              'rlnAmplitudeContrast',
                              'rlnAngleRot',
                              'rlnAngleTilt',
                              'rlnAnglePsi',
                              'rlnCtfBfactor',
                              'rlnCtfScalefactor',
                              ])

    @ staticmethod
    def _getImportChoices():
        """ Return a list of possible choices
        to choose to which kind of volume the CTF3D should be estimated for.
        """
        return ['Per volume', 'Per subvolume']

    def _estimateCTF3DPerSubvolume(self, tsExp, setTsInfo, tsCounter):
        starFileList = []
        mrcFileList = []
        tiltList = []
        doseList = []
        coordCounter = 0
        ts = tsExp.getTS()
        ctfs = tsExp.getCTFSeries()
        coordList = tsExp.getCoords()
        if not self._doseFromMdoc:
            tiltList = tsExp.getTiltAngles()
            doseList = tsExp.getDoses()

        sizeX, _, sizeZ, _ = ImageHandler().getDimensions(tsExp.getTS().getFirstItem().getFileName())
        for coord in coordList:
            tomoTable = self._createTable()
            ctf3DStar = self._getCtfFile(ts.getTsId(),
                                         coordCounter=coordCounter,
                                         fileExt=CTFSTAR,
                                         ctfMode=CTF3D_PER_SUBVOLUME)
            mrc3DStar = self._getCtfFile(ts.getTsId(),
                                         coordCounter=coordCounter,
                                         fileExt=CTFMRC,
                                         ctfMode=CTF3D_PER_SUBVOLUME)
            coordCounter += 1
            starFileList.append(ctf3DStar)
            mrcFileList.append(mrc3DStar)
            self.ctfMRCFileList.append(mrc3DStar)

            for ti, ctf in zip(ts, ctfs):
                avgDefocus = (ctf.getDefocusU() + ctf.getDefocusV()) / 2
                tiltAngleDegs = ti.getTiltAngle()
                tiltAngleRads = np.deg2rad(tiltAngleDegs)
                xTomo = float(coord.getX() - (sizeX / 2)) * setTsInfo[SRATE]
                zTomo = float(coord.getZ() - (sizeZ / 2)) * setTsInfo[SRATE]
                # Calculating the height difference of the particle from the tilt axis
                xImg = (xTomo * (math.cos(tiltAngleRads))) + (zTomo * (math.sin(tiltAngleRads)))
                deltaD = xImg * math.sin(tiltAngleRads)
                partDef = avgDefocus + deltaD
                # Weighting the 3D CTF model using the tilt dependent scale factor and the dose dependent B-Factor
                tiltScale = math.cos(abs(tiltAngleRads))
                if self._doseFromMdoc:
                    tiltImgDose = ti.getAcquisition().getDosePerFrame()
                else:
                    tiltImgDose = self._getCurrentDose(tiltAngleDegs, tiltList, doseList)
                doseWeight = tiltImgDose * self.bFactor
                # Add row to table
                tomoTable.addRow(partDef,
                                 setTsInfo[VOLTAGE],
                                 setTsInfo[SPHAB],
                                 setTsInfo[AMPCN],
                                 0.0,
                                 tiltAngleDegs,
                                 0.0,
                                 doseWeight,
                                 tiltScale)
            # Write the corresponding CTF star file
            tomoTable.write(ctf3DStar)

        self.tsExpandedList[tsCounter].setCTFStarList(starFileList)
        self.tsExpandedList[tsCounter].setCTFMRCList(mrcFileList)

    def _estimateCTF3DPerVolume(self, tsExp, setTsInfo, tsCounter):
        starFileList = []
        mrcFileList = []
        tiltList = []
        doseList = []
        ts = tsExp.getTS()
        ctfs = tsExp.getCTFSeries()
        if not self._doseFromMdoc:
            tiltList = tsExp.getTiltAngles()
            doseList = tsExp.getDoses()

        tomoTable = self._createTable()
        ctf3DStar = self._getCtfFile(ts.getTsId(),
                                     fileExt=CTFSTAR,
                                     ctfMode=CTF3D_PER_VOLUME)
        mrc3DStar = self._getCtfFile(ts.getTsId(),
                                     fileExt=CTFMRC,
                                     ctfMode=CTF3D_PER_VOLUME)
        starFileList.append(ctf3DStar)
        mrcFileList.append(mrc3DStar)
        self.ctfMRCFileList.append(mrc3DStar)

        for ti, ctf in zip(ts, ctfs):
            avgDefocus = (ctf.getDefocusU() + ctf.getDefocusV()) / 2
            tiltAngleDegs = ti.getTiltAngle()
            tiltAngleRads = np.deg2rad(tiltAngleDegs)
            # Weighting the 3D CTF model using the tilt dependent scale factor and the dose dependent B-Factor
            tiltScale = math.cos(abs(tiltAngleRads))
            if self._doseFromMdoc:
                tiltImgDose = ti.getAcquisition().getDosePerFrame()
            else:
                tiltImgDose = self._getCurrentDose(tiltAngleDegs, tiltList, doseList)
            doseWeight = tiltImgDose * self.bFactor
            # Add row to table
            tomoTable.addRow(avgDefocus,
                             setTsInfo[VOLTAGE],
                             setTsInfo[SPHAB],
                             setTsInfo[AMPCN],
                             0.0,
                             tiltAngleDegs,
                             0.0,
                             doseWeight,
                             tiltScale)

        tomoTable.write(ctf3DStar)
        self.tsExpandedList[tsCounter].setCTFStarList(starFileList)
        self.tsExpandedList[tsCounter].setCTFMRCList(mrcFileList)


class ExtendedTS:
    """This class represents a expanded version of the tilt series, and adds the additional data
    required to calculate the CTF3d"""
    def __init__(self, ts, ctfs, doseFileObj, coords):
        self._ts = ts
        self._ctfs = ctfs
        self._doseFile = doseFileObj
        self._coords = coords
        self._ctfStarFileList = []
        self._ctfMrcFileList = []

    def getTS(self):
        return self._ts

    def getCTFSeries(self):
        return self._ctfs

    def getDoseFile(self):
        return self._doseFile

    def getTiltAngles(self):
        if self.getDoseFile():
            return self.getDoseFile().getTiltAngles()
        else:
            return [ti.getTiltAngle() for ti in self.getTS()]

    def getDoses(self):
        if self.getDoseFile():
            return self.getDoseFile().getDoses()
        else:
            return [ti.getAcquisition().getDosePerFrame() for ti in self.getTS()]

    def getCoords(self):
        return self._coords

    def getCTFStarList(self):
        return self._ctfStarFileList

    def setCTFStarList(self, ctfStarList):
        self._ctfStarFileList = ctfStarList

    def getCTFMRCList(self):
        return self._ctfMrcFileList

    def setCTFMRCList(self, mrcStarList):
        self._ctfMrcFileList = mrcStarList


# This class represents a dose file, and provides getters for its contents
class DoseFile:
    def __init__(self, doseFile):
        self._fileName = doseFile
        self._tiltAngles = []
        self._doses = []

        self._readDoseFile(doseFile)

    def _readDoseFile(self, doseFile):
        fid = open(doseFile, 'r')
        lines = fid.readlines()
        fid.close()
        angles = []
        doses = []
        for line in lines:
            lineCols = line.split()
            angles.append(float(lineCols[0]))
            doses.append(float(lineCols[1]))

        self._tiltAngles = angles
        self._doses = doses

    def getTiltAngles(self):
        return self._tiltAngles

    def getDoses(self):
        return self._doses

    def getDoseFileNmae(self):
        return self._fileName
