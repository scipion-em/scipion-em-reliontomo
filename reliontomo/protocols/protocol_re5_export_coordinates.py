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
import enum
import logging
from typing import Union
import numpy as np
from pwem.convert import transformations
from pwem.protocols import EMProtocol
from pyworkflow import BETA
from pyworkflow.object import Pointer
from pyworkflow.protocol import PointerParam, IntParam, LEVEL_ADVANCED
from pyworkflow.utils import Message, cyanStr, redStr
from reliontomo.convert.convertBase import getTransformInfoFromCoordOrSubtomo
from reliontomo.objects import RelionSetOfPseudoSubtomograms, RelionPSubtomogram
from tomo.constants import SCIPION
from tomo.objects import SetOfCoordinates3D, SetOfTomograms, Coordinate3D

logger = logging.getLogger(__name__)
IN_PARTICLES = 'inputSubTomos'
IN_TOMOS = 'inTomos'

class ExportCoordsRe5(enum.Enum):
    coordinates3d = SetOfCoordinates3D


class ProtTomoExportRe5Coords(EMProtocol):
    """
    Extracts and exports 3D particle coordinates from Relion pseudo-
    subtomograms into tomogram-based coordinate sets. The protocol is
    designed to recover particle positions and orientations so they can
    be reused in downstream tomography workflows, visualization, or
    re-extraction procedures.

    AI Generated:

    Export 3D Coordinates from Relion Particles
    (ProtTomoExportRe5Coords) - User Manual

        Overview

        The Export 3D Coordinates from Relion Particles protocol converts
        Relion pseudo-subtomogram information into a standard set of 3D
        coordinates associated with tomograms. Its primary purpose is to
        recover the spatial localization of particles after subtomogram
        extraction, classification, refinement, or other tomography
        processing stages.

        In cryo-electron tomography workflows, this protocol is
        particularly useful when users need to return from processed
        particles back to their original biological context inside the
        tomogram. This allows the recovered coordinates to be reused for
        re-extraction with different box sizes, visualization of particle
        distributions, structural interpretation inside cellular
        environments, or preparation of downstream subtomogram averaging
        workflows.

        Biological Context and Typical Applications

        In many tomography projects, particles undergo several rounds of
        classification, refinement, and cleaning before biologically
        relevant subsets are identified. Once meaningful particle
        populations have been selected, researchers often need to map
        those particles back into the tomogram in order to inspect their
        native spatial organization.

        This protocol enables that transition from processed particle
        datasets back to tomographic coordinates. Typical applications
        include visualization of ribosomes, membrane complexes, viral
        particles, cytoskeletal assemblies, or macromolecular machines
        directly within reconstructed cellular volumes.

        Another common application involves re-extraction of particles
        using larger box sizes, different binning factors, or updated
        reconstruction parameters. This is especially important when
        downstream refinement requires higher-resolution information or
        when earlier extraction settings were optimized only for rapid
        exploratory analysis.

        Inputs and Dataset Compatibility

        The protocol requires two principal inputs: a set of Relion
        pseudo-subtomograms and a corresponding set of tomograms. These
        datasets must share compatible tomogram identifiers so that each
        particle can be correctly associated with its biological source
        volume.

        Dataset consistency is biologically important because incorrect
        associations between particles and tomograms may produce
        misplaced coordinates, incorrect orientations, or biologically
        meaningless spatial interpretations.

        In practice, users should verify that the tomograms correspond
        to the same reconstruction stage and sampling conditions used
        during particle extraction and refinement.

        Coordinate Scaling and Sampling Considerations

        The protocol automatically adapts particle coordinates to the
        sampling properties of the target tomograms. This step is
        essential because pseudo-subtomograms and tomograms may exist at
        different pixel sizes or binning levels.

        Correct coordinate scaling ensures that exported particles are
        accurately positioned within the reconstructed cellular volume.
        From a biological perspective, proper scaling is critical for
        meaningful interpretation of spatial organization, molecular
        crowding, membrane association, or intracellular localization.

        The protocol also allows specification of the coordinate box
        size. This parameter defines the approximate particle region
        associated with each exported coordinate and becomes especially
        important during later re-extraction or visualization workflows.

        Orientation Recovery and Spatial Context

        In addition to particle positions, the protocol preserves
        orientation information associated with the pseudo-subtomograms.
        This allows exported coordinates to retain their rotational
        context relative to the tomogram.

        Preserving orientation information is biologically valuable when
        studying ordered assemblies, membrane-associated complexes,
        filamentous systems, or anisotropic macromolecular distributions.
        It also facilitates downstream alignment and refinement workflows
        that depend on maintaining consistent particle geometry.

        Outputs and Interpretation

        The protocol produces a set of 3D coordinates directly linked to
        the selected tomograms. These coordinates can be visualized in
        tomographic viewers, reused for extraction procedures, or passed
        into additional cryo-electron tomography processing pipelines.

        Each exported coordinate preserves both positional and
        transformation information, allowing particles to remain linked
        to their refined structural interpretation while being expressed
        in the tomogram reference frame.

        Biologically, the resulting coordinate set allows users to
        investigate how particles are spatially organized within the
        specimen. This can reveal clustering behavior, membrane
        localization, intracellular gradients, lattice organization, or
        structural heterogeneity across different cellular regions.

        Practical Recommendations

        Before running the protocol, users should confirm that the
        tomograms and particle datasets correspond to the same biological
        specimen and reconstruction workflow. Mismatched datasets may
        generate incomplete or incorrect coordinate assignments.

        When planning re-extraction workflows, it is generally advisable
        to verify that the selected box size is large enough to contain
        the full particle and surrounding structural context, especially
        for flexible or membrane-associated complexes.

        For visualization purposes, inspecting the exported coordinates
        directly within tomographic viewers is strongly recommended. This
        helps validate that particles are correctly positioned and that
        biologically meaningful spatial patterns are preserved.

        Final Perspective

        Exporting refined particle information back into tomographic
        space is an important step in many cryo-electron tomography
        studies because it reconnects high-quality structural results
        with their native biological environment. By preserving particle
        localization and orientation within tomograms, this protocol
        enables more complete interpretation of molecular organization,
        cellular architecture, and in situ structural biology.
    """

    _label = 'export 3D coordinates from Relion particles'
    _devStatus = BETA
    _possibleOutputs = ExportCoordsRe5

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.tomosDict = None
        
    # --------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        form.addSection(label=Message.LABEL_INPUT)

        form.addParam(IN_PARTICLES, PointerParam,
                      pointerClass=RelionSetOfPseudoSubtomograms,
                      label='Relion particles',
                      important=True)

        form.addParam(IN_TOMOS, PointerParam,
                      pointerClass=SetOfTomograms,
                      label='Tomograms',
                      important=True,
                      help='Select the tomograms to which you want to\n'
                           'associate the coordinates.')

        form.addParam('boxSize', IntParam,
                      allowsNull=True,
                      expertLevel=LEVEL_ADVANCED,
                      label='Box Size (opt)',
                      help='Determine the box size of the extracted coordinates. By default, '
                           'the program assigns the box size directly from the coordinates '
                           'associated to the subtomograms.')

    # --------------------------- INSERT steps functions ----------------------
    def _insertAllSteps(self):
        self._initialize()
        self._insertFunctionStep(self.exportCoordsStep, needsGPU=False)

    def _initialize(self):
        inTomograms = self.getInputTomograms()
        inParticles = self.getInputParticles()

        # Compute matching TS id among coordinates, the tilt-series and the CTFs, they all could be a subset
        particlesTsIds = set(inParticles.getTSIds())
        tomoTsIds = set(inTomograms.getTSIds())
        presentTsIds = particlesTsIds & tomoTsIds
        nonMatchingTsIds = (particlesTsIds ^ tomoTsIds) - presentTsIds
        # Validate the intersection
        if len(presentTsIds) <= 0:
            raise Exception("There isn't any common tsIds among the tomograms and the "
                            "particles introduced.")
        if len(nonMatchingTsIds) > 0:
            logger.info(cyanStr(f"TsIds not common in the introduced tomograms and "
                                f"particles are: {nonMatchingTsIds}"))

        self.tomosDict = {tomo.getTsId(): tomo.clone() for tomo in inTomograms.iterItems()
                          if tomo.getTsId() in presentTsIds}

    def exportCoordsStep(self):
        inParticlesPointer = self.getInputParticles(returnPointer=True)
        inTomogramsPointer = self.getInputTomograms(returnPointer=True)
        inParticles = inParticlesPointer.get()
        partSRate = inParticles.getSamplingRate()
        scaleFactor = self._getCoordsScaleFactor()
        boxSize = int(inParticles.getBoxSize() * scaleFactor)
        boxSize = boxSize if boxSize % 2 == 0 else boxSize + 1
        tomosSRate = inTomogramsPointer.get().getSamplingRate()
        outAttribName = self._possibleOutputs.coordinates3d.name

        outCoords = SetOfCoordinates3D.create(self.getPath(), template='coordinates%s.sqlite')
        outCoords.setSamplingRate(tomosSRate)
        outCoords.setPrecedents(self.getInputTomograms(returnPointer=True))
        outCoords.setBoxSize(boxSize)

        for tsId, tomo in self.tomosDict.items():
            logger.info(cyanStr(f'tsId = {tsId} - tomogram {tomo}. Exporting the coordinates...'))
            for pSubtomo in inParticles.iterSubtomos(volume=tomo):
                try:
                    coord = Coordinate3D()
                    coord.setVolume(tomo)
                    x = scaleFactor * pSubtomo.getXInImg() / partSRate  # Originally in angstrom at the scale of the particle
                    y = scaleFactor * pSubtomo.getYInImg() / partSRate
                    z = scaleFactor * pSubtomo.getZInImg() / partSRate
                    coord.setPosition(x, y, z, SCIPION)  # Already centered in Relion 5
                    trMatrix = self.getTransformMatrix(pSubtomo, sRate=partSRate, scaleFactor=scaleFactor)
                    coord.setMatrix(trMatrix)
                    coord.setGroupId(pSubtomo.getGroupId())
                    outCoords.append(coord)

                except Exception as e:
                    logger.error(redStr(f'tsId = {tsId} - psubtomo {pSubtomo} failed: {e}. Skipping...'))
                    continue

        if len(outCoords) == 0:
            raise Exception(f'No output/s {outAttribName} were generated. '
                            f'Please check the Output Log > run.stdout and run.stderr')

        self._defineOutputs(**{outAttribName: outCoords})
        self._defineSourceRelation(inTomogramsPointer, outCoords)
        self._defineSourceRelation(inParticlesPointer, outCoords)

    # --------------------------- UTILS functions -----------------------------
    def getInputParticles(self, returnPointer: bool=False) -> Union[Pointer, RelionSetOfPseudoSubtomograms]:
        reParticlesPointer = getattr(self, IN_PARTICLES)
        return reParticlesPointer if returnPointer else reParticlesPointer.get()

    def getInputTomograms(self, returnPointer: bool=False) -> Union[Pointer, SetOfTomograms]:
        inTomosPointer = getattr(self, IN_TOMOS)
        return inTomosPointer if returnPointer else inTomosPointer.get()

    @staticmethod
    def getTransformMatrix(particle: RelionPSubtomogram,
                           sRate: float = 1.,
                           scaleFactor: float = 1.) -> np.ndarray:

        # angles[0],  # 2, rlnTomoSubtomogramRot
        # angles[1],  # 3, rlnTomoSubtomogramTilt
        # angles[2],  # 4, rlnTomoSubtomogramPsi
        # pSubtomo.getRot(),  # 5, rlnAngleRot
        # pSubtomo.getTilt(),  # 6, rlnAngleTilt
        # pSubtomo.getPsi(),  # 7, rlnAnglePsi
        # pSubtomo.getTiltPrior(),  # 8, rlnAngleTiltPrior
        # pSubtomo.getPsiPrior(),  # 9, rlnAnglePsiPrior
        # shifts[0],  # 14, rlnOriginXAngst
        # shifts[1],  # 15, rlnOriginYAngst
        # shifts[2],  # 16, rlnOriginZAngst,
        anglesTomoSubtomo, _ = getTransformInfoFromCoordOrSubtomo(particle, sRate)
        shiftx = scaleFactor * particle.getX() / sRate  # They are stored in Relion 5 in angstroms, but needed here in px
        shifty = scaleFactor * particle.getY() / sRate
        shiftz = scaleFactor * particle.getZ() / sRate
        rot = particle.getRot()
        tilt = particle.getTilt()
        psi = particle.getPsi()
        angles = (rot, tilt, psi)
        shifts = np.array([shiftx, shifty, shiftz])

        radAnglesTomoSubtomo = -np.deg2rad(anglesTomoSubtomo)
        radAngles = -np.deg2rad(angles)
        MPicking = transformations.euler_matrix(*radAnglesTomoSubtomo, 'szyz')
        MRefine = transformations.euler_matrix(*radAngles, 'szyz')
        # MRefine[0, 3] = -shifts[0]
        # MRefine[1, 3] = -shifts[1]
        # MRefine[2, 3] = -shifts[2]
        M = np.dot(MPicking, MRefine)
        M = np.linalg.inv(M)

        # rotShifts = - np.linalg.inv(MPicking[:3,:3]).dot(shifts)
        # rotShifts = - MPicking[:3,:3].dot(shifts)
        # M[0, 3] = rotShifts[0]
        # M[1, 3] = rotShifts[1]
        # M[2, 3] = rotShifts[2]
        return M

    def _getCoordsScaleFactor(self) -> float:
        return self.getInputParticles().getSamplingRate() / self.getInputTomograms().getSamplingRate()

    # -------------------------- INFO functions -------------------------------
    def _summary(self):
        msg = []
        if self.isFinished():
            scaleFactor = self._getCoordsScaleFactor()
            msg.append(f'The particles introduced were scaled using an scale factor of '
                       f'{scaleFactor:.2f} to be expressed considering '
                       f'the size of the introduced tomograms')

        return msg
