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
from os import mkdir
from os.path import exists
from pyworkflow.utils import getParentFolder
from reliontomo.constants import SUBTOMO_NAME, FILE_NOT_FOUND, PSUBTOMOS_SQLITE
from reliontomo.convert import createReaderTomo
from reliontomo.objects import RelionSetOfPseudoSubtomograms
from reliontomo.protocols.protocol_base_import_from_star import ProtBaseImportFromStar
from reliontomo.protocols.protocol_base_relion import IS_RELION_50
from reliontomo.utils import getAbsPath
from tomo.objects import SetOfSubTomograms, SetOfCoordinates3D

SUBTOMOGRAMS_DIR = 'Subtomograms'

class outputObjects(Enum):
    coordinates = SetOfCoordinates3D
    subtomograms = SetOfSubTomograms


class ProtImportSubtomogramsFromStar(ProtBaseImportFromStar):
    """
    Imports subtomogram datasets from STAR metadata files and integrates them
    into tomography processing workflows. The protocol is intended to recover
    previously generated subtomograms or pseudo-subtomograms together with their
    associated tomographic context, allowing users to continue analysis within
    a unified cryo-electron tomography environment.

    AI Generated:

    Import Subtomograms From STAR (ProtImportSubtomogramsFromStar) - User Manual
        Overview

        The Import Subtomograms From STAR protocol is designed to incorporate
        subtomogram datasets described in STAR metadata tables into an existing
        tomography project. In cryo-electron tomography workflows, STAR files
        are commonly used to store references to reconstructed particles,
        acquisition information, and associated coordinate metadata. This
        protocol provides a structured way to reconnect those metadata
        descriptions with the tomographic context required for downstream
        analysis.

        The protocol is particularly useful when subtomograms have been
        generated outside the current project, transferred from another RELION
        workflow, or archived independently from the original processing
        session. By importing the STAR-based metadata, users can restore the
        relationship between reconstructed particles and their originating
        tomograms without manually rebuilding the dataset structure.

        Inputs and Biological Context

        The protocol requires a STAR file containing references to subtomograms
        or pseudo-subtomograms, together with a corresponding set of source
        tomograms. The tomograms provide the acquisition geometry, sampling
        information, and spatial reference system necessary to interpret the
        imported particles correctly.

        From a biological perspective, preserving the connection between
        subtomograms and their parent tomograms is essential. Many downstream
        analyses, including subtomogram averaging, classification, particle
        localization, and structural interpretation, depend on accurate spatial
        context. Improperly linked metadata may lead to orientation errors,
        incorrect scaling, or loss of biological meaning.

        Compatibility With Different RELION Generations

        The protocol supports both traditional subtomogram representations and
        newer pseudo-subtomogram workflows introduced in modern RELION
        tomography pipelines. This flexibility allows users to integrate data
        originating from different processing generations within a consistent
        project structure.

        Traditional subtomograms generally represent reconstructed 3D particle
        volumes extracted from tomograms, whereas pseudo-subtomograms are
        intermediate representations optimized for RELION tomography refinement
        strategies. Supporting both formats is particularly important in
        facilities or collaborative projects where datasets may originate from
        different software environments or processing standards.

        Validation and Data Integrity

        One of the most important aspects of the protocol is verification of
        the metadata and file structure referenced by the STAR file. The import
        process ensures that the required fields are present and that the files
        referenced in the metadata can be located correctly.

        In practical cryo-ET workflows, metadata inconsistencies are a common
        source of failure. STAR files may reference moved files, renamed
        directories, or incomplete exports from external systems. Detecting
        these issues early prevents downstream processing errors that may be
        difficult to diagnose later.

        Biological users should pay particular attention to voxel size
        consistency and tomogram provenance. Even when subtomograms appear to
        import correctly, differences in sampling rate or acquisition geometry
        between projects can affect averaging accuracy and structural
        interpretation.

        Outputs and Their Interpretation

        After completion, the protocol generates a structured subtomogram or
        pseudo-subtomogram dataset linked to the corresponding tomographic
        information. The imported particles preserve their original metadata
        associations and become immediately available for subsequent processing
        stages such as classification, refinement, alignment, or visualization.

        The resulting dataset serves as a bridge between external reconstruction
        workflows and integrated tomography analysis environments. This is
        particularly valuable in collaborative cryo-ET studies where particle
        extraction and refinement may occur across multiple computational
        platforms.

        Practical Recommendations

        Before import, users should verify that all files referenced in the
        STAR metadata remain accessible and that directory structures have not
        changed since export. Maintaining consistent project organization
        greatly reduces import problems and ensures reproducible workflows.

        When combining datasets originating from different microscopes,
        acquisition sessions, or processing pipelines, users should carefully
        confirm that voxel sizes and coordinate conventions remain compatible.
        Small inconsistencies at this stage can propagate into significant
        biological interpretation errors during averaging or classification.

        For large collaborative projects, it is often advisable to preserve the
        original STAR files together with their associated particle directories
        as archival records. This simplifies later re-import and improves
        reproducibility across institutions or software versions.

        Final Perspective

        Importing subtomogram metadata is more than a technical data-loading
        operation. It is a critical step in preserving the biological and
        spatial relationships that underpin cryo-electron tomography analysis.
        Careful validation of metadata consistency, sampling information, and
        tomographic provenance ensures that imported datasets remain reliable
        for structural interpretation and downstream refinement workflows.
    """

    _label = 'import subtomograms from a star file'
    _possibleOutputs = outputObjects

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.isCoordsFile = False

    # --------------------------- STEPS functions -----------------------------

    def _initialize(self):
        super()._initialize()
        # Generate the firectoy in which the linked subtomograms pointed from the star file will be stored
        mkdir(self._getExtraPath(SUBTOMOGRAMS_DIR))

    def _importStep(self):
        precedentsSet = self.inTomos.get()
        # Generate the corresponding precedents and 3d coordinates
        super()._importStep()
        # Generate the set of subtomograms
        acq = precedentsSet.getAcquisition()
        sRate = precedentsSet.getSamplingRate()
        if self.readerVersion < 5:
            subtomoSet = SetOfSubTomograms.create(self._getPath(), template='setOfSubTomograms%s.sqlite')
            subtomoSet.setSamplingRate(sRate)
            subtomoSet.setAcquisition(acq)
            self.reader.starFile2SubtomogramsImport(subtomoSet,
                                                    getattr(self, outputObjects.coordinates.name),
                                                    self._getExtraPath(SUBTOMOGRAMS_DIR),
                                                    getParentFolder(self.starFile.get()))
        else:
            subtomoSet = RelionSetOfPseudoSubtomograms.create(self.getPath(), template=PSUBTOMOS_SQLITE)
            subtomoSet.setSamplingRate(sRate)
            subtomoSet.setAcquisition(acq)
            self.reader.starFile2PseudoSubtomograms(subtomoSet)
        self._defineOutputs(**{outputObjects.subtomograms.name: subtomoSet})
        self._defineSourceRelation(self.inTomos, subtomoSet)

    # --------------------------- INFO functions -----------------------------
    def _validate(self):
        errorMsg = super()._validate()
        reader, isReader40 = createReaderTomo(self.starFile.get(), isRelion5=False)
        errorsInPointedFiles = self._checkFilesPointedFromStarFile(getParentFolder(self.starFile.get()),
                                                                   reader.dataTable)
        if errorsInPointedFiles:
            errorMsg.append(errorsInPointedFiles)
        return errorMsg

    # --------------------------- UTILS functions ------------------------------
    def _checkFilesPointedFromStarFile(self, starFilePath, dataTable):
        errorsFound = ''
        fields2check = [SUBTOMO_NAME]
        filesPattern = '\tRow %i - %s\n'

        errorsFound += self._checkFieldsInDataTable(dataTable, fields2check)
        # Check if the files pointed from those fields exist
        if errorsFound:
            filesErrorMsgHead = 'The following files were not found [row, tomoFile, subtomoFile]:\n'
            for counter, row in enumerate(dataTable):
                subtomoFileNotFound = self._fileNotFound(row, SUBTOMO_NAME, starFilePath)
                if subtomoFileNotFound:
                    errorsFound += filesPattern % (counter, subtomoFileNotFound)

            if errorsFound:
                errorsFound = filesErrorMsgHead + errorsFound

        return errorsFound

    @staticmethod
    def _checkFieldsInDataTable(dataTable, fieldList):
        fieldErrors = ''
        if fieldList:
            fieldNotFoundPattern = 'Fields %s were not found in the star file introduced.\n'
            notFoundFields = [field for field in fieldList if not dataTable.hasColumn(field)]
            if notFoundFields:
                pattern = '[%s]' % (' '.join(notFoundFields))
                fieldErrors = (fieldNotFoundPattern % pattern)

        return fieldErrors

    @staticmethod
    def _fileNotFound(row, field, starFilePath):
        statusMsg = ''
        tomoFile = row.get(field, FILE_NOT_FOUND)
        tomoFileAbs = getAbsPath(starFilePath, tomoFile)
        if not exists(tomoFileAbs):
            statusMsg = tomoFile

        return statusMsg
