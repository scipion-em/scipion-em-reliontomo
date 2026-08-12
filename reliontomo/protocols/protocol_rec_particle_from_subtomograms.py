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
from pwem.convert.headers import fixVolume
from pyworkflow.protocol.params import (PointerParam, FloatParam, StringParam, EnumParam, IntParam, LEVEL_ADVANCED)
from pwem.protocols import ProtReconstruct3D
from reliontomo.protocols.protocol_base_relion import ProtRelionTomoBase, IS_RELION_50
from tomo.objects import AverageSubTomogram
from reliontomo.convert import writeSetOfSubtomograms


class outputObjects(Enum):
    average = AverageSubTomogram


class ProtRelionSubTomoReconstructAvg(ProtReconstruct3D):
    """
    Reconstructs an averaged three-dimensional volume from a set of
    aligned subtomograms using Relion reconstruction strategies. The
    protocol combines particle information from multiple tomographic
    observations to generate a consensus structural map suitable for
    visualization, interpretation, and downstream cryo-electron
    tomography analysis.

    AI Generated:

    Average Subtomogram Reconstruction (ProtRelionSubTomoReconstructAvg)
    — User Manual

        Overview

        This protocol reconstructs a three-dimensional averaged volume from
        a collection of aligned subtomograms using Relion-based
        reconstruction methods. Its primary purpose is to generate a
        consensus structural representation by combining information from
        multiple particles extracted from cryo-electron tomography datasets.

        In biological cryo-electron tomography workflows, subtomogram
        averaging is one of the most important strategies for improving the
        signal-to-noise ratio of macromolecular structures embedded in
        native environments. By averaging many aligned particles, the
        protocol enhances reproducible structural features while reducing
        random experimental noise.

        The reconstructed volume can subsequently be used for biological
        interpretation, visualization, refinement, classification, or
        integration with complementary structural biology workflows. Typical
        applications include the study of membrane complexes, viral
        assemblies, ribosomes, cytoskeletal filaments, and molecular
        machines within cellular environments.

        Inputs and General Workflow

        The protocol requires a set of aligned subtomograms as input. These
        subtomograms should already contain reliable orientation and
        positional information obtained from previous alignment or
        refinement procedures. The reconstruction quality depends strongly
        on the consistency and accuracy of these alignments.

        During reconstruction, the protocol combines projection information
        from all selected particles into a common three-dimensional
        coordinate system. The resulting map represents the average
        structural state of the analyzed particle population.

        In practical biological workflows, reconstruction quality improves
        as the number of accurately aligned particles increases. However,
        strong structural heterogeneity within the dataset may reduce map
        sharpness and obscure biologically meaningful conformational states.

        Symmetry Considerations

        The protocol allows the application of molecular symmetry during
        reconstruction. Correct symmetry selection can substantially improve
        map quality because equivalent structural views are averaged
        together to enhance the effective signal-to-noise ratio.

        Symmetry is especially beneficial for highly ordered biological
        assemblies such as viral capsids, symmetric oligomers, or molecular
        cages. Applying the proper symmetry may reveal structural details
        that remain difficult to observe in asymmetric reconstructions.

        However, incorrect symmetry assignment can introduce severe
        artifacts and distort biologically relevant asymmetry or
        conformational variability. For heterogeneous or asymmetric
        complexes, reconstruction without symmetry constraints is generally
        recommended.

        The protocol follows Relion symmetry conventions, allowing users to
        define common point-group symmetries appropriate for the biological
        system under investigation. More info:
        http://www2.mrc-lmb.cam.ac.uk/relion/index.php/Conventions_%26_File_formats#Symmetry

        Resolution Control and Fourier Space Processing

        The protocol allows control over the maximum spatial resolution
        considered during reconstruction. Restricting the reconstruction to
        lower spatial frequencies may improve robustness in noisy datasets
        or early exploratory analyses.

        In biological practice, conservative resolution limits are often
        useful during initial reconstruction stages, particularly when the
        alignment accuracy remains uncertain. Higher-resolution refinement
        should generally be attempted only after stable alignment has been
        confirmed.

        Fourier-space padding can also be adjusted to improve interpolation
        accuracy during reconstruction. Increased padding generally improves
        numerical precision and reduces interpolation artifacts, although it
        requires additional computational resources.

        In many routine workflows, default padding values provide a good
        balance between computational efficiency and reconstruction quality.

        Half-Set Reconstruction and Validation

        The protocol supports reconstruction of independent particle
        subsets, including half-set reconstructions commonly used for
        validation and resolution estimation. These independent
        reconstructions are biologically important because they allow
        objective assessment of structural reproducibility.

        Half-map comparisons are frequently used in Fourier Shell
        Correlation analysis to estimate the effective resolution of the
        reconstruction. Consistent agreement between independently refined
        half-maps indicates that the recovered structural information is
        reproducible rather than dominated by noise.

        The ability to reconstruct specific subsets also facilitates focused
        biological analyses, including comparison of structural states,
        validation of classification strategies, or investigation of
        dataset heterogeneity.

        Class-Based Reconstruction

        The protocol can optionally restrict reconstruction to particles
        belonging to a selected structural class. This functionality is
        biologically valuable when heterogeneous datasets contain multiple
        conformational states, assembly intermediates, or distinct molecular
        populations.

        Focused reconstruction of a single class often improves structural
        homogeneity and reveals biologically meaningful conformational
        differences that may otherwise be obscured in global averages.

        Careful class selection is particularly important when studying
        flexible assemblies, dynamic molecular machines, or complexes with
        variable ligand occupancy.

        Computational Considerations

        Three-dimensional subtomogram reconstruction is computationally
        intensive due to the large volume of tomographic data and the
        extensive Fourier-space calculations required for averaging. The
        protocol supports parallel processing strategies to improve
        computational efficiency.

        Computational cost depends strongly on factors such as particle
        number, box size, reconstruction resolution, and applied symmetry.
        Larger datasets and higher-resolution targets require substantially
        greater processing time and memory resources.

        In practical workflows, it is often advantageous to begin with
        lower-resolution exploratory reconstructions before proceeding to
        more demanding high-resolution analyses.

        Outputs and Biological Interpretation

        The protocol produces an averaged three-dimensional subtomogram map
        suitable for visualization and downstream analysis. The resulting
        reconstruction represents the consensus structural information
        recovered from the selected particle population.

        Improved averaging often reveals biologically important features
        such as secondary structure organization, domain architecture,
        membrane association, ligand density, or macromolecular assembly
        organization.

        The reconstructed map may subsequently serve as input for further
        refinement, classification, segmentation, docking, or comparative
        structural analysis workflows.

        Biological interpretation should always consider the degree of
        structural heterogeneity present in the input dataset. Averaging
        heterogeneous conformations may blur flexible regions and obscure
        meaningful dynamic behavior.

        Practical Recommendations

        In most biological workflows, the quality of subtomogram alignment
        is the single most important factor determining reconstruction
        quality. Users should carefully validate alignment consistency
        before attempting high-resolution averaging.

        Symmetry should only be imposed when supported by independent
        biological evidence. Incorrect symmetry assignment may produce maps
        that appear visually sharp but lack biological validity.

        When analyzing heterogeneous datasets, reconstructing separate
        particle classes often yields more interpretable biological results
        than combining all particles into a single average.

        Conservative resolution limits and moderate padding settings are
        generally appropriate during early exploratory analyses. More
        aggressive reconstruction settings should be introduced gradually as
        confidence in the dataset quality increases.

        Final Perspective

        Subtomogram averaging reconstruction is a foundational technique in
        cryo-electron tomography because it enables recovery of structural
        information from extremely noisy experimental data acquired within
        native biological environments.

        Careful alignment, biologically appropriate symmetry selection,
        objective validation through half-map analysis, and thoughtful
        treatment of structural heterogeneity are essential for generating
        reliable and biologically meaningful three-dimensional
        reconstructions.
    """
    _label = 'Average subtomo'
    inStarName = 'input_particles'
    outTomoName = 'output_volume'
    _possibleOutputs = outputObjects

    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        form.addSection(label='Input')

        form.addParam('inputSubtomos', PointerParam,
                      pointerClass='SetOfSubTomograms',
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
        form.addParallelSection(threads=0, mpi=1)

    # -------------------------- INSERT steps functions -----------------------
    def _insertAllSteps(self):
        self._initialize()
        self._insertFunctionStep(self.convertInputStep)
        self._insertFunctionStep(self.reconstructStep)
        self._insertFunctionStep(self.createOutputStep)

    def _getProgram(self, program='relion_reconstruct'):
        """ Get the program name depending on the MPI use or not. """
        if self.numberOfMpi > 1:
            program += '_mpi'
        return program

    # -------------------------- STEPS functions ------------------------------
    def _initialize(self):
        """ Centralize how files are called for iterations and references. """
        myDict = {
            self.inStarName: self._getExtraPath(self.inStarName + '.star'),
            self.outTomoName: self._getExtraPath(self.outTomoName + '.mrc')
            }
        self._updateFilenamesDict(myDict)

    def convertInputStep(self):
        """ Create the input file in STAR format as expected by Relion."""
        subtomosSet = self.inputSubtomos.get()
        starFile = self._getFileName(self.inStarName)
        # This binary is expecting the star file column names in relion 3 format
        writeSetOfSubtomograms(subtomosSet, starFile)

    def reconstructStep(self):
        self.runJob(self._getProgram(), self._genReconstructCmd())

    def createOutputStep(self):
        imgSet = self.inputSubtomos.get()
        volume = AverageSubTomogram()
        volumeFile = self._getFileName(self.outTomoName)
        fixVolume(volumeFile)
        volume.setFileName(volumeFile)
        volume.setSamplingRate(imgSet.getSamplingRate())
        
        self._defineOutputs(**{outputObjects.average.name: volume})
        self._defineSourceRelation(self.inputSubtomos, volume)
    
    # -------------------------- INFO functions -------------------------------
    def _validate(self):
        pass
    
    def _summary(self):
        summary = []
        if not hasattr(self, 'outputVolume'):
            summary.append("Output tomogram not ready yet.")
        else:
            summary.append("Output tomogram has been reconstructed.")

        return summary

    # --------------------------- UTILS functions -----------------------------
    @classmethod
    def isDisabled(cls):
        """ Return True if this Protocol is disabled.
        Disabled protocols will not be offered in the available protocols."""
        return True if IS_RELION_50 else False

    def _genReconstructCmd(self):
        imgSet = self.inputSubtomos.get()

        params = ' --i %s' % self._getFileName(self.inStarName)
        params += ' --o %s' % self._getFileName(self.outTomoName)
        params += ' --sym %s' % self.symmetryGroup.get()
        params += ' --angpix %0.3f' % imgSet.getSamplingRate()
        params += ' --maxres %0.3f' % self.maxRes.get()
        params += ' --pad %0.3f' % self.pad.get()

        subset = -1 if self.subset.get() == 0 else self.subset
        params += ' --subset %d' % subset

        # if Plugin.IS_GT30():
        params += ' --class %d' % self.classNum.get()

        if self.extraParams.hasValue():
            params += " " + self.extraParams.get()

        return params
