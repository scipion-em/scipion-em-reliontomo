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

from pyworkflow import BETA, PROD
from reliontomo import Plugin
from reliontomo.protocols.protocol_base_import_from_star import ProtBaseImportFromStar
from tomo.objects import SetOfTomograms, SetOfCoordinates3D

IS_RELION_50 = Plugin.isRe50()


class outputObjects(Enum):
    tomograms = SetOfTomograms()
    coordinates = SetOfCoordinates3D()


class ProtImportCoordinates3DFromStar(ProtBaseImportFromStar):
    """
    Imports three-dimensional particle coordinates from RELION STAR metadata
    files and integrates them into cryo-electron tomography workflows. The
    protocol is intended to recover spatial particle locations identified in
    external processing pipelines and make them available for downstream
    subtomogram extraction, refinement, classification, and structural analysis.

    AI Generated:

    Import 3D Coordinates From STAR (ProtImportCoordinates3DFromStar) - User Manual
        Overview

        The Import 3D Coordinates From STAR protocol provides a mechanism for
        incorporating particle coordinate information stored in RELION STAR
        files into a tomography project. In cryo-electron tomography workflows,
        particle coordinates define the spatial positions of biological
        structures within reconstructed tomograms and serve as the foundation
        for subtomogram extraction and averaging procedures.

        This protocol is especially useful when coordinates have been generated
        externally, transferred from another project, or produced by automated
        particle-picking strategies outside the current workflow environment.
        By importing the coordinate metadata directly, users can continue
        analysis without manually recreating particle positions.

        Biological Context and Coordinate Interpretation

        In cryo-ET, coordinates correspond to the estimated locations of
        macromolecular complexes, ribosomes, membrane proteins, viral particles,
        or other biological assemblies embedded within the cellular or
        vitrified environment. Accurate coordinate definition is essential
        because every downstream operation depends on the spatial precision of
        these particle locations.

        Even small coordinate inaccuracies can affect subtomogram centering,
        alignment stability, and the final biological interpretation of
        averaged structures. For this reason, imported coordinates should be
        interpreted within the context of the original tomographic acquisition,
        voxel size calibration, and particle-picking strategy.

        Inputs and Workflow Integration

        The protocol requires a RELION-compatible STAR file containing 3D
        coordinate information. These coordinates are associated with source
        tomograms that define the acquisition geometry and spatial reference
        frame. Maintaining this relationship is critical for ensuring that the
        imported particles remain biologically meaningful and geometrically
        consistent.

        Imported coordinates can subsequently be used for subtomogram
        extraction, pseudo-subtomogram generation, alignment workflows, or
        classification pipelines. This makes the protocol an important bridge
        between particle detection and structural interpretation.

        Compatibility and Experimental Status

        The protocol supports coordinate import workflows associated with modern
        RELION tomography pipelines and is particularly relevant for projects
        integrating RELION 5 tomography metadata. Because tomography metadata
        conventions continue to evolve, users should carefully verify imported
        coordinate consistency before large-scale downstream processing.

        In collaborative environments, differences in software versions,
        coordinate conventions, or tomogram orientation standards may influence
        how imported coordinates are interpreted. Careful validation is
        therefore recommended when combining data generated by different
        laboratories or processing facilities.

        Outputs and Their Interpretation

        After execution, the protocol produces a structured set of 3D
        coordinates associated with their corresponding tomograms and metadata.
        These imported coordinates become part of the tomography project and
        can be directly used by subsequent subtomogram analysis protocols.

        The resulting coordinate dataset preserves the spatial organization of
        the biological particles inside the tomographic volume. This enables
        downstream analyses that depend on accurate localization, including
        structural averaging, spatial distribution studies, and contextual
        interpretation within cellular environments.

        Practical Recommendations

        Before importing coordinates, users should verify that the tomograms
        referenced in the project correspond exactly to those used during the
        original particle-picking process. Differences in voxel size, binning,
        cropping, or tomogram orientation can introduce systematic coordinate
        offsets that compromise downstream refinement.

        It is also advisable to inspect a subset of imported coordinates
        visually within the tomographic volume before proceeding to extensive
        processing. Early validation helps identify scaling mismatches,
        coordinate inversion problems, or metadata inconsistencies that might
        otherwise propagate through the workflow.

        For large cryo-ET studies, maintaining consistent coordinate conventions
        across all processing stages significantly improves reproducibility and
        reduces the risk of biological misinterpretation.

        Final Perspective

        Importing 3D coordinates is a foundational step in tomography analysis
        because it establishes the spatial relationship between detected
        particles and the biological specimen. Reliable coordinate integration
        ensures that subsequent subtomogram extraction and averaging procedures
        remain geometrically accurate and biologically interpretable throughout
        the cryo-electron tomography workflow.
    """

    _label = 'import 3D coordinates from a star file'
    _devStatus = BETA if IS_RELION_50 else PROD
    _possibleOutputs = outputObjects

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.isCoordsFile = True
