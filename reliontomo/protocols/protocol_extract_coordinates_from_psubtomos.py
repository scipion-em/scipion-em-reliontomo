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
from pyworkflow import BETA
from pyworkflow.object import Integer
from reliontomo.constants import RELION_3D_COORD_ORIGIN
from reliontomo.convert.convertBase import genTransformMatrix
from reliontomo.protocols.protocol_base_refine import ProtRelionRefineBase
from tomo.objects import SetOfTomograms, SetOfCoordinates3D, Coordinate3D
from tomo.utils import getObjFromRelation


class outputObjects(Enum):
    coordinates = SetOfCoordinates3D()


class ProtExtractCoordsFromPSubtomos(ProtRelionRefineBase):
    """Protocol to extract a set of 3D coordinates from a set of pseudo-subtomograms"""

    _label = 'extract coordinates from pseudo-subtomograms'
    _possibleOutputs = outputObjects
    _devStatus = BETA

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.tomoSet = None

    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        super()._defineCommonInputParams(form)

    # -------------------------- INSERT steps functions -----------------------
    def _insertAllSteps(self):
        self._initialize()
        self._insertFunctionStep(self._extractCoordsStep)

    # -------------------------- STEPS functions ------------------------------
    def _initialize(self):
        # Get the tomograms through the relations
        self.tomoSet = getObjFromRelation(self.inReParticles.get(), self, SetOfTomograms)

    def _extractCoordsStep(self):
        pSubtomos = self.inReParticles.get()
        # Manage the precedents
        precedents = self.tomoSet
        precedentsSRate = precedents.getSamplingRate()
        precedentIdDict = {}
        for tomo in precedents:
            precedentIdDict[tomo.getTsId()] = tomo.clone()

        # Generate the output set of coordinates
        coordSet = SetOfCoordinates3D.create(self._getPath(), template='coordinates%s.sqlite')
        coordSet.setPrecedents(precedents)
        coordSet.setBoxSize(pSubtomos.getBoxSize())
        coordSet.setSamplingRate(precedentsSRate)
        scaleFactor = self.getScaleFactor()

        for pSubtomo in pSubtomos:
            # CHeck if the current psubtomo has a precedent in the tomoset
            coordTsId = pSubtomo.getTsId()
            precedent = precedentIdDict.get(coordTsId, None)
            if precedent:
                coordinate3d = Coordinate3D()
                coordinate3d.setTomoId(coordTsId)
                coordinate3d.setVolume(precedent)
                # Scale x, y and z properly
                # NOTE: pSubtomos and coordinate3D are both Scipion EMObjects, so no conversion between spaces
                # is required here
                coordinate3d.setX(pSubtomo.getX() * scaleFactor, RELION_3D_COORD_ORIGIN)
                coordinate3d.setY(pSubtomo.getY() * scaleFactor, RELION_3D_COORD_ORIGIN)
                coordinate3d.setZ(pSubtomo.getZ() * scaleFactor, RELION_3D_COORD_ORIGIN)
                # Scale the shifts properly
                tm = pSubtomo.getTransform().getMatrix()
                tm[0, 3] = tm[0, 3] * scaleFactor
                tm[1, 3] = tm[0, 1] * scaleFactor
                tm[2, 3] = tm[0, 2] * scaleFactor
                coordinate3d.setMatrix(self.getTransformMatrix(pSubtomo, precedentsSRate))
                # Extended fields
                coordinate3d._classNumber = Integer(pSubtomo.getClassId())
                coordinate3d._randomSubset = Integer(pSubtomo.getRdnSubset())
                # Append the current coordinate to the output set
                coordSet.append(coordinate3d)

        self._defineOutputs(**{outputObjects.coordinates.name: coordSet})
        self._defineSourceRelation(pSubtomos, coordSet)

    # -------------------------- INFO functions -------------------------------

    def _summary(self):
        summary = []
        if self.isFinished():
            outputCoords = getattr(self, outputObjects.coordinates.name, None)
            if outputCoords:
                summary.append('Coordinates were scaled to the sampling rate of the tomograms from where\n'
                               'they were picked --> *%.2f Å/pix*.' % outputCoords.getSamplingRate())

        return summary

    # --------------------------- UTILS functions -----------------------------

    def getScaleFactor(self):
        inPSubtomos = self.inReParticles.get()
        sRatePSubtomo = inPSubtomos.getSamplingRate() / inPSubtomos.getRelionBinning()
        sRateTomo = self.tomoSet.getSamplingRate()
        return sRatePSubtomo / sRateTomo

    @staticmethod
    def getTransformMatrix(pSubtomo, tomoSRate, invert=True):
        shiftx = pSubtomo._sxAngst.get() / tomoSRate
        shifty = pSubtomo._syAngst.get() / tomoSRate
        shiftz = pSubtomo._szAngst.get() / tomoSRate
        tilt = pSubtomo._rot.get()
        psi = pSubtomo._tilt.get()
        rot = pSubtomo._psi.get()

        return genTransformMatrix(shiftx, shifty, shiftz, rot, tilt, psi, invert)
