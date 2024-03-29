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
from os.path import exists, join
from emtable import Table

from pwem.objects import VolumeMask, FSC
from pwem.protocols import EMProtocol
from pyworkflow.protocol import PointerParam, StringParam
from pyworkflow.utils import Message, createLink
from reliontomo.constants import IN_PARTICLES_STAR, POSTPROCESS_DIR, OPTIMISATION_SET_STAR, PSUBTOMOS_SQLITE, \
    OUT_PARTICLES_STAR
from reliontomo.convert import writeSetOfPseudoSubtomograms, readSetOfPseudoSubtomograms
from reliontomo.objects import createSetOfRelionPSubtomograms, RelionSetOfPseudoSubtomograms


class ProtRelionTomoBase(EMProtocol):

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    @staticmethod
    def _defineCommonInputParams(form):
        form.addSection(label=Message.LABEL_INPUT)
        form.addParam('inReParticles', PointerParam,
                      important=True,
                      pointerClass='RelionSetOfPseudoSubtomograms',
                      label='Pseudo-Subtomograms',
                      help='Pseudo-subtomograms do not aim to accurately represent the scattering potential of '
                           'the underlying particles. Instead, they serve as a practical means to implement an '
                           ' approximation to the 2D approach within the existing RELION framework. In the original'
                           ' RELION4 article an accurate defition is given, see: \n '
                           'https://doi.org/10.7554/eLife.83724 \n '
                           ' A more technical explanation, pseudo-subtomograms are 3D-Arrays (volumes) that constructed '
                           ' from the sums of 2D tilt-series images pre-multiplied by contrast transfer functions (CTFs), '
                           ' along with auxiliary arrays that store the corresponding sum of squared CTFs and the frequency '
                           ' of observation for each 3D voxel.')

    @staticmethod
    def _defineExtraParams(form, addAdditionalSection=True):
        if addAdditionalSection:
            form.addSection(label='Additional')
        form.addParam('extraParams', StringParam,
                      label='Additional arguments',
                      help="In this box command-line arguments may be provided that are not generated by the GUI. This "
                           "may be useful for testing developmental options and/or expert use of the program, e.g: \n"
                           "--verb 1\n"
                           "--pad 2\n")

    def getInputParticles(self):
        return self.inReParticles.get()

    def getOutStarFileName(self):
        return self._getExtraPath(IN_PARTICLES_STAR)

    def genInStarFile(self, withPriors=False):
        """It will check if the set size and the stored particles star file are of the same size or not. In
        the first case, a link will be made to the previous particles star file to avoid generating it and in the
        second case, a new file will be generated containing only the ones present in the input set."""
        inReParticlesSet = self.inReParticles.get()
        outStarFileName = self.getOutStarFileName()
        if inReParticlesSet.getSize() == inReParticlesSet.getNReParticles() and not withPriors:
             self.info("Using existing star (%s) file instead of generating a new one." % inReParticlesSet.getParticles())
             createLink(inReParticlesSet.getParticles(), outStarFileName)
        else:
            writeSetOfPseudoSubtomograms(inReParticlesSet, outStarFileName, withPriors=withPriors)

    def _genPostProcessOutputMrcFile(self, fileName):
        """File generated using the sharpening protocol (called post-process protocol) and also using the
        rec particle from TS protocol in case the optional input 'solvent mask' is introduced."""
        postProccesMrc = VolumeMask()
        postProccesMrc.setFileName(self._getExtraPath(POSTPROCESS_DIR, fileName))
        postProccesMrc.setSamplingRate(self.inReParticles.get().getCurrentSamplingRate())

        return postProccesMrc

    def genRelionParticles(self, optimisationFileName=OPTIMISATION_SET_STAR, particles=OUT_PARTICLES_STAR, binningFactor=None, boxSize=24):
        """Generate a RelionSetOfPseudoSubtomograms object containing the files involved for the next protocol,
        considering that some protocols don't generate the optimisation_set.star file. In that case, the input Object
        which represents it will be copied and, after that, this method will be used to update the corresponding
        attribute."""

        # Create the set
        inParticlesSet = self.inReParticles.get()
        psubtomoSet = RelionSetOfPseudoSubtomograms.create(self.getPath(), template=PSUBTOMOS_SQLITE)
        psubtomoSet.copyInfo(inParticlesSet)

        # Verify out star file
        extraPath = self._getExtraPath()
        optimSetStar = join(extraPath, optimisationFileName)
        if exists(optimSetStar):
            psubtomoSet.filesMaster = optimSetStar

        particles = join(extraPath, particles)
        if exists(particles):
            psubtomoSet._particles.set(particles)


        if binningFactor:
            psubtomoSet.setRelionBinning(binningFactor)

        if boxSize:
            psubtomoSet.setBoxSize(boxSize)

        # Fill the items (pseudo subtomos/particles) from de particles star file
        psubtomoSet.setSamplingRate(psubtomoSet.getCurrentSamplingRate())
        readSetOfPseudoSubtomograms(psubtomoSet)

        return psubtomoSet

    def genFSCs(self, starFile, tableName, fscColumns):
        fscSet = self._createSetOfFSCs()
        table = Table(fileName=starFile, tableName=tableName)
        resolution_inv = table.getColumnValues('rlnResolution')
        for columnName in fscColumns:
            columValues = table.getColumnValues(columnName)
            fsc = FSC(objLabel=columnName[3:])
            fsc.setData(resolution_inv, columValues)
            fscSet.append(fsc)

        fscSet.write()
        return fscSet

    def _genExtraParamsCmd(self):
        return ' ' + self.extraParams.get() if self.extraParams.get() else ''

