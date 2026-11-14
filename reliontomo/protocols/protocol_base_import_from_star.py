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
from enum import Enum
from os.path import exists
from emtable import Table
from pwem.protocols import EMProtocol
from pyworkflow.object import Boolean
from pyworkflow.protocol import FileParam, FloatParam, IntParam, PointerParam
from pyworkflow.utils import Message, getParentFolder, createLink, makePath
from reliontomo import Plugin
from reliontomo.constants import TOMO_NAME, IN_PARTICLES_STAR, PIXEL_SIZE, TOMO_PARTICLE_ID
from reliontomo.convert import createReaderTomo
from reliontomo.protocols.protocol_base_relion import IS_RELION_50
from tomo.protocols.protocol_base import ProtTomoBase
from tomo.objects import SetOfCoordinates3D

logger = logging.getLogger(__file__)
IS_RE5_PICKING_ATTR = '_relion5Picking'


class importCoordsOutputs(Enum):
    coordinates = SetOfCoordinates3D


class ProtBaseImportFromStar(EMProtocol, ProtTomoBase):
    """
    Imports tomographic particle coordinates from STAR metadata files and associates them with an existing set of tomograms.

    AI Generated:

    Import Coordinates From STAR (ProtBaseImportFromStar) - User Manual
        Overview

        The Import Coordinates From STAR protocol is designed to integrate particle coordinate information generated in
        external cryo-electron tomography workflows into a Scipion tomography project. Its primary purpose is to read
        coordinates stored in STAR metadata files and connect them with previously imported tomograms so they can be
        used in downstream subtomogram analysis, particle extraction, refinement, or classification workflows.

        In practical cryo-ET processing, coordinates are often generated outside the main processing environment using
        particle picking tools, template matching approaches, neural-network detectors, or manual annotation software.
        This protocol provides a bridge between those coordinate sources and the tomography data already available in
        the project. The resulting coordinate set becomes a structured object that can be reused consistently across
        multiple downstream protocols.

        Inputs and Coordinate Association

        The protocol requires two principal inputs: a STAR metadata file containing coordinate information and a set of
        tomograms that define the spatial reference system. The imported coordinates are interpreted relative to the
        geometry and sampling characteristics of the selected tomograms. Correct correspondence between the STAR file
        identifiers and the tomogram identifiers is therefore essential for successful import.

        From a biological perspective, this association step ensures that each detected particle is linked to the
        correct cellular volume or tilt series reconstruction. Incorrect matching between coordinates and tomograms may
        lead to particles being misplaced in unrelated regions, producing biologically meaningless downstream analyses.

        Sampling Rate Considerations

        One of the most important practical aspects of coordinate import is the management of sampling rates. The
        protocol allows coordinates generated at one pixel size to be correctly scaled into tomograms reconstructed at
        another pixel size. This situation is common in cryo-ET workflows where particle picking is performed on
        downsampled tomograms for speed, while later subtomogram processing may use higher-resolution reconstructions.

        If a coordinate sampling rate is not explicitly provided, the protocol attempts to infer it from the metadata
        available in the STAR file. When no such information exists, the tomogram sampling rate is assumed. Correct
        scaling is biologically critical because inaccurate coordinate scaling may shift particles away from their true
        structural locations, compromising extraction quality and subsequent refinement.

        Compatibility Across STAR Formats

        Cryo-ET metadata formats have evolved across different versions of RELION and related tomography tools. This
        protocol is intended to support coordinate import from multiple generations of STAR files while maintaining
        compatibility with modern tomography workflows. Older coordinate conventions and newer RELION tomography data
        structures are both considered during import validation.

        In practical use, this flexibility is particularly valuable for laboratories combining historical datasets with
        newly acquired experiments or integrating coordinates produced by independent external software. The protocol
        attempts to identify the metadata conventions automatically and adapt the import process accordingly.

        Coordinate Validation and Tomogram Matching

        Before importing coordinates, the protocol performs consistency checks between the coordinate metadata and the
        tomogram collection. These validations help ensure that identifiers used in the STAR file correspond to the
        tomograms available in the project. This step is especially important in large cryo-ET studies where multiple
        tomograms or experimental conditions are processed simultaneously.

        Biologically, proper validation reduces the risk of cross-assignment between unrelated tomograms, which could
        otherwise produce misleading localization patterns or incorrect subtomogram populations. Careful identifier
        consistency is therefore essential for reliable interpretation.

        Box Size and Extraction Context

        The protocol allows the user to define a coordinate box size that represents the expected extraction region
        surrounding each particle. This parameter influences the spatial context retained around the imported particle
        coordinates and therefore affects later subtomogram reconstruction and averaging stages.

        Smaller box sizes may be appropriate for compact macromolecular assemblies, while larger values are often
        beneficial when studying membrane-associated complexes, flexible structures, or crowded intracellular
        environments where contextual information may assist alignment and classification.

        Outputs and Their Interpretation

        After execution, the protocol produces a set of 3D coordinates linked directly to the corresponding tomograms.
        These coordinates preserve the spatial organization of the original tomography experiment and can be used in
        particle extraction, subtomogram averaging, classification, or visualization workflows.

        The resulting coordinate set serves as a foundational object within cryo-ET processing pipelines because many
        downstream analyses depend on accurate particle localization. The biological quality of later reconstructions is
        therefore strongly influenced by the precision and consistency of the imported coordinates.

        Practical Recommendations

        In routine practice, users should verify that tomogram identifiers in the STAR file match the identifiers used
        in the tomography dataset before import. It is also advisable to confirm that coordinate scaling corresponds to
        the intended reconstruction level, particularly when coordinates were generated on binned tomograms.

        When integrating coordinates from older software environments, careful visual inspection after import is highly
        recommended to ensure that particles occupy biologically plausible regions within the tomograms. Small metadata
        inconsistencies can otherwise propagate into substantial downstream alignment errors.

        Final Perspective

        Coordinate import is a critical integration step in cryo-electron tomography because it connects external
        particle detection workflows with the internal processing environment used for subtomogram analysis. Accurate
        association between coordinates, tomograms, and sampling rates ensures that downstream structural refinement
        begins from biologically meaningful particle locations and preserves the integrity of the original experiment.
    """

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.linkedStarFileName = IN_PARTICLES_STAR
        self.linkedTomosDirName = 'Tomograms'
        self.reader = None
        self.readerVersion = None
        self.coordsSRate = None
        self.starFilePath = None
        self.isCoordsFile = None  # To be defined by the child classes

    @classmethod
    def isDisabled(cls):
        """ Return True if this Protocol is disabled.
        Disabled protocols will not be offered in the available protocols."""
        return True if IS_RELION_50 else False

    def _defineParams(self, form):
        form.addSection(label=Message.LABEL_INPUT)
        form.addParam('starFile', FileParam,
                      label='Star file')
        form.addParam('inTomos', PointerParam,
                      pointerClass='SetOfTomograms',
                      label='Tomograms',
                      help='Tomograms to which the coordinates will be referred.')
        form.addParam('samplingRate', FloatParam,
                      label='Coordinates sampling rate [Å/pix] (opt.)',
                      allowsNull=True,
                      help='If empty, it will be considered to be the same as the tomograms. The ratio of both '
                           'tomograms and coordinates sampling rate will be used to scale the coordinates '
                           'properly to the tomograms introduced if they are different.\n\n'
                           'For OLD versions of the star file:'
                           'If empty, the protocol will try to read it from the label "_r%s" '
                           'if it is present in the introduced star file. If not find, it will be '
                           'considered to be the same as the tomograms. The ratio of both tomograms '
                           'and coordinates sampling rate will be used to scale the coordinates '
                           'properly to the tomograms introduced if they are different' % PIXEL_SIZE)
        form.addParam('boxSize', IntParam,
                      label='Box Size [pix]',
                      default=20)

    # --------------------------- STEPS functions -----------------------------

    def _insertAllSteps(self):
        self._initialize()
        self._insertFunctionStep(self._importStep, needsGPU=False)

    def _initialize(self):
        # Get the star file and make a symbolic link in extra directory
        self.starFilePath = getParentFolder(self.starFile.get())
        newStarName = self._getExtraPath(self.linkedStarFileName)
        createLink(self.starFile.get(), newStarName)
        # Read the star file
        self.reader, self.readerVersion = createReaderTomo(starFile=newStarName,
                                                           isRelion5=IS_RELION_50,
                                                           isCoordsStar=self.isCoordsFile)
        # Generate the directory in which the linked tomograms pointed from the star file will be stored
        makePath(self._getExtraPath(self.linkedTomosDirName))
        # Get the coordinates sampling rate
        self.coordsSRate = self.getCoordsSRate()

    def _importStep(self):
        inTomosPointer = self.inTomos
        precedentsSet = inTomosPointer.get()
        inTomosSRate = precedentsSet.getSamplingRate()
        coordSet = SetOfCoordinates3D.create(self._getPath(), template='coordinates%s.sqlite')
        coordSet.setPrecedents(inTomosPointer)  # As a pointer is better for streaming
        coordSet.setSamplingRate(inTomosSRate)  # They will be re-scaled to the tomograms sRate
        coordSet.setBoxSize(self.boxSize.get())
        self.reader.starFile2Coords3D(coordSet, precedentsSet, self.coordsSRate / inTomosSRate)
        setattr(coordSet, IS_RE5_PICKING_ATTR, Boolean(self.checkIfRe5PickedCoords()))

        self._defineOutputs(**{importCoordsOutputs.coordinates.name: coordSet})
        self._defineSourceRelation(inTomosPointer, coordSet)

    # --------------------------- INFO functions ------------------------------

    def _validate(self):
        errors = []
        # Check if the files referred in the star file exists
        starFile = self.starFile.get()
        if not exists(starFile):
            errors.append('It was not possible to locate the introduced file. Please check the path.\n%s' % starFile)
            return errors
        # Check the compatibility between the introduced file and the version of Relion used by the plugin
        isRe5Star = self.checkIfRe5PickedCoords()
        reader, readerVersion = createReaderTomo(self.starFile.get(), isRelion5=IS_RELION_50)
        if self.isCoordsFile:
            if isRe5Star:
                if Plugin.isRe40():
                    errors.append('The introduced star file is in Relion 5 format, while the plugin is configured to '
                                  'use Relion 4.')
                else:
                    self._checkIfTsIdErrorsInStar(reader, errors)
            elif readerVersion == 4:
                self._checkIfTsIdErrorsInStar(reader, errors)
            # Older star files like the ones generated by pyseg cannot be checked like this because they do not
            # contain the tsId
        return errors

    def _summary(self):
        summary = []
        if self.isFinished():
            if not self.samplingRate.get():
                summary.append('The sampling rate considered was the one read from the star file.')
            summary.append("Coordinates imported from %s." % self.starFile.get())

        return summary

    # --------------------------- UTILS functions ------------------------------
    def checkIfRe5PickedCoords(self) -> bool:
        dataTable = Table()
        dataTable.read(self.starFile.get())
        presentColumns = dataTable.getColumnNames()
        # The column psSegImage is present in PySeg star files (picked coordinates that follow Relion 3 format)
        # The column rlnTomoParticleId is present in Relion 4 data model, but not in Relion 5
        if 'psSegImage' in presentColumns or TOMO_PARTICLE_ID in presentColumns:
            return False
        else:
            return True
        # """
        # In Pelion 5 file particles.star, there's a table named "general" that contains one field and one values to
        # indicate if the particles are 2D or 3D. Value 1 is used for 2D particles while value 0 is used for 3D
        # particles. Example:
        # data_general
        # _rlnTomoSubTomosAre2DStacks                       0
        # """
        # try:
        #     reader, _ = createReaderTomo(self.starFile.get(), tableName=GENERAL_TABLE)
        #     return bool(reader.dataTable.getColumnValues(RLN_ARE2DSTACKS)[0])
        # except Exception:
        #     return False  # That table does not exist in star files generated with Relion versions older than 5

    def getCoordsSRate(self):
        if self.samplingRate.get():
            coordsSRate = self.samplingRate.get()
        else:
            sRateFromStar = self.reader.dataTable[0].get(PIXEL_SIZE, None)
            if sRateFromStar:
                coordsSRate = float(sRateFromStar)
            else:
                coordsSRate = self.inTomos.get().getSamplingRate()
        return coordsSRate

    def _checkIfTsIdErrorsInStar(self, reader, errorMsg: list):
        tsIds = [tomo.getTsId() for tomo in self.inTomos.get()]
        if tsIds:
            coordTomoIds = list(set([row.get(TOMO_NAME) for row in reader.dataTable]))
            matchingCoorTomoIds = [coordTomoId for coordTomoId in coordTomoIds if coordTomoId in tsIds]
            if len(set(coordTomoIds) & set(matchingCoorTomoIds)) == 0:
                errorMsg.append(f'No matching TsIds were found between the\n'
                                f'\t- Coordinates rlnTomoName {coordTomoIds}\n'
                                f'\t- Tilt series tsId {matchingCoorTomoIds}]')
        else:
            errorMsg.append('TsId is empty in the introduced tilt series. It is not possible to match the coordinates '
                            'and the tilt series in that case.')
