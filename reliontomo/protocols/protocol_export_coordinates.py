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
from pyworkflow.object import Pointer, Float
from pyworkflow.protocol import PointerParam, IntParam, LEVEL_ADVANCED
from pyworkflow.utils import Message, cyanStr, redStr
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
    Extract the coordinates information from a set of pseudo-subtomograms.

    This protocol is useful when we want to re-extract the subtomograms
    (maybe resulting from classification) with the
    original dimensions. It can be also handy to visualize the resulting
    subtomograms in their location on the tomograms.
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

        # The coordinates need to be re-scaled to be at the same size of the tomograms
        self.coordsScaleFactor = Float(inParticles.getSamplingRate() / inTomograms.getSamplingRate())
        self._store(self.coordsScaleFactor)

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
        tomosSRate = inTomogramsPointer.get().getSamplingRate()
        outAttribName = self._possibleOutputs.coordinates3d.name
        outCoords = SetOfCoordinates3D.create(self.getPath(), template='coordinates%s.sqlite')
        self._defineOutputs(**{outAttribName: outCoords})
        self._defineSourceRelation(inTomogramsPointer, outCoords)
        self._defineSourceRelation(inParticlesPointer, outCoords)
        outCoords.setSamplingRate(tomosSRate)
        outCoords.setPrecedents(self.getInputTomograms(returnPointer=True))
        partSRate = inParticles.getSamplingRate()
        scaleFactor = self.coordsScaleFactor.get()

        for tsId, tomo in self.tomosDict.items():
            logger.info(cyanStr(f'tsId = {tsId} - tomogram {tomo}. Exporting the coordinates...'))
            for pSubtomo in inParticles.iterSubtomos(volume=tomo):
                try:
                    logger.info(cyanStr(f'--> {pSubtomo}'))
                    coord = Coordinate3D()
                    coord.setVolume(tomo)
                    x = scaleFactor * pSubtomo.getXInImg() / partSRate  # Originally in agstrom at the scale of the particle
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
        shiftx = scaleFactor * particle.getX() / sRate  # They are stored in Relion 5 in angstroms, but needed here in px
        shifty = scaleFactor * particle.getY() / sRate
        shiftz = scaleFactor * particle.getZ() / sRate
        rot = particle.getRot()
        tilt = particle.getTilt()
        psi = particle.getPsi()
        angles = (rot, tilt, psi)
        shifts = (shiftx, shifty, shiftz)
        radAngles = -np.deg2rad(angles)
        M = transformations.euler_matrix(radAngles[0], radAngles[1], radAngles[2], 'szyz')
        M[0, 3] = -shifts[0]
        M[1, 3] = -shifts[1]
        M[2, 3] = -shifts[2]
        M = np.linalg.inv(M)
        return M

    # -------------------------- INFO functions -------------------------------
    def _summary(self):
        msg = []
        if self.isFinished():
            msg.append(f'The particles introduced were scaled using an scale factor of '
                       f'{self.coordsScaleFactor.get():.2f} to be expressed considering '
                       f'the size of the introduced tomograms')

        return msg
