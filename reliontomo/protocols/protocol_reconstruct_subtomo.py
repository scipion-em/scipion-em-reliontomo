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

from pyworkflow.protocol.params import (PointerParam, FloatParam,  
                                        StringParam, BooleanParam,
                                        EnumParam, IntParam, LEVEL_ADVANCED)

from pwem.protocols import ProtReconstruct3D
from tomo.objects import Tomogram, AverageSubTomogram
from relion import Plugin

from reliontomo.convert import writeSetOfSubtomograms


class ProtRelionSubTomoReconstruct(ProtReconstruct3D):
    """ This protocol reconstructs a volume using Relion.

    Reconstruct a volume from a given set of particles.
    The alignment parameters will be converted to a Relion star file
    and used as direction projections to reconstruct.
    """
    _label = 'tomo reconstruct'
    inStarName = 'input_particles'
    outTomoName = 'output_volume'

    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        form.addSection(label='Input')

        form.addParam('inputSubtomos', PointerParam,
                      pointerClass='SetOfSubTomograms',
                      # pointerCondition='hasAlignmentProj',
                      label="Input subtomograms",
                      help='Select the input subtomograms from the project.')
        form.addParam('symmetryGroup', StringParam, default='c1',
                      label="Symmetry group",
                      help='See [[Relion Symmetry][http://www2.mrc-lmb.cam.ac.uk/'
                           'relion/index.php/Conventions_%26_File_formats#Symmetry]] '
                           'page for a description of the symmetry format '
                           'accepted by Relion')
        form.addParam('maxRes', FloatParam, default=-1,
                      label="Maximum resolution (A)",  
                      help='Maximum resolution (in Angstrom) to consider \n'
                           'in Fourier space (default Nyquist).')
        form.addParam('pad', FloatParam, default=2,
                      label="Padding factor")
        form.addParam('subset', EnumParam, default=0,
                      choices=['all', 'half1', 'half2'],
                      display=EnumParam.DISPLAY_HLIST,
                      label='Subset to reconstruct',
                      help='Subset of images to consider.')
        # if Plugin.IS_GT30():
        form.addParam('classNum', IntParam, default=-1,
                      label='Use only this class',
                      help='Consider only this class (-1: use all classes)')
        
        form.addParam('extraParams', StringParam, default='',
                      expertLevel=LEVEL_ADVANCED,
                      label='Extra parameters: ', 
                      help='Extra parameters to *relion_reconstruct* program. '
                           'Address to Relion to see full list of options.')
        form.addSection('CTF')
        form.addParam('doCTF', BooleanParam, default=False,
                      label='Apply CTF correction?')
        form.addParam('ctfIntactFirstPeak', BooleanParam, default=False,
                      condition='doCTF',
                      label='Leave CTFs intact until first peak?')
        
        form.addParallelSection(threads=0, mpi=1)

    # -------------------------- INSERT steps functions -----------------------
    def _insertAllSteps(self):
        self._createFilenameTemplates()
        self._insertFunctionStep('convertInputStep')
        self._insertReconstructStep()
        self._insertFunctionStep('createOutputStep')

    def _getProgram(self, program='relion_reconstruct'):
        """ Get the program name depending on the MPI use or not. """
        if self.numberOfMpi > 1:
            program += '_mpi'
        return program

    def _insertReconstructStep(self):
        imgSet = self.inputSubtomos.get()

        params = ' --i %s' % self._getFileName(self.inStarName)
        params += ' --o %s' % self._getFileName(self.outTomoName)
        params += ' --sym %s' % self.symmetryGroup.get()
        params += ' --angpix %0.3f' % imgSet.getSamplingRate()
        params += ' --maxres %0.3f' % self.maxRes.get()
        params += ' --pad %0.3f' % self.pad.get()

        subset = -1 if self.subset.get() == 0 else self.subset
        params += ' --subset %d' % subset

        if Plugin.IS_GT30():
            params += ' --class %d' % self.classNum.get()

        if self.doCTF:
            params += ' --ctf'
            if self.ctfIntactFirstPeak:
                params += ' --ctf_intact_first_peak'

        if self.extraParams.hasValue():
            params += " " + self.extraParams.get()

        self._insertFunctionStep('reconstructStep', params)

    # -------------------------- STEPS functions ------------------------------
    def reconstructStep(self, params):
        self.runJob(self._getProgram(), params)

    def _createFilenameTemplates(self):
        """ Centralize how files are called for iterations and references. """
        myDict = {
            'input_particles': self._getTmpPath(self.inStarName + '.star'),
            'output_volume': self._getExtraPath(self.outTomoName + '.mrc')
            }
        self._updateFilenamesDict(myDict)

    def convertInputStep(self):
        """ Create the input file in STAR format as expected by Relion.
        """
        # TODO: If the input particles comes from Relion, just link the file.
        imgSet = self.inputSubtomos.get()
        imgStar = self._getFileName(self.inStarName)
        writeSetOfSubtomograms(imgSet, imgStar, outputDir=self._getTmpPath())

    def createOutputStep(self):
        imgSet = self.inputSubtomos.get()
        volume = AverageSubTomogram()
        volume.setFileName(self._getFileName(self.outTomoName))
        volume.setSamplingRate(imgSet.getSamplingRate())
        
        self._defineOutputs(outputAvgSubtomogram=volume)
        self._defineSourceRelation(self.inputSubtomos, volume)
    
    # -------------------------- INFO functions -------------------------------
    def _validate(self):
        errors = []

        return errors
    
    def _summary(self):
        summary = []
        if not hasattr(self, 'outputVolume'):
            summary.append("Output tomogram not ready yet.")
        else:
            summary.append("Output tomogram has been reconstructed.")

        return summary
