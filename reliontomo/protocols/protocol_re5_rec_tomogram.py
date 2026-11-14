# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors:     Scipion Team
# *
# * National Center of Biotechnology, CSIC, Spain
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
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************
from enum import Enum

from pyworkflow.protocol import PointerParam, IntParam, FloatParam, BooleanParam, EnumParam, StringParam
from pyworkflow.utils import Message
from reliontomo import Plugin
from reliontomo.constants import IN_TS_STAR, TOMOGRAMS_DIR
from reliontomo.convert import convert50_tomo
from reliontomo.protocols.protocol_base_relion import ProtRelionTomoBase, IS_RELION_50
from reliontomo.utils import getProgram
from tomo.objects import SetOfTomograms, Tomogram, TiltSeries

# Reconstruct options
SINGLE_TOMO = 0
ALL_TOMOS = 1


class tomoRecOtputObjects(Enum):
    tomograms = SetOfTomograms


class ProtRelion5TomoReconstruct(ProtRelionTomoBase):
    """
    Reconstructs cryo-electron tomography volumes from aligned tilt
    series and their associated CTF information. The protocol generates
    tomograms suitable for visualization, particle picking, subtomogram
    analysis, denoising workflows, and downstream structural studies.

    AI Generated:

    Tomograms Reconstruction (ProtRelion5TomoReconstruct) - User Manual

        Overview

        The Tomograms Reconstruction protocol reconstructs 3D tomograms
        from aligned tilt-series datasets using Relion tomography tools.
        Its primary purpose is to transform a series of 2D projection
        images acquired at multiple tilt angles into biologically
        interpretable 3D cellular or molecular volumes. These tomograms
        provide the structural context necessary for subtomogram
        averaging, particle extraction, in situ structural biology, and
        visual inspection of macromolecular organization inside cells or
        vitrified specimens.

        In practical cryo-electron tomography workflows, this protocol
        is often one of the central reconstruction stages following
        motion correction, tilt-series alignment, and CTF estimation.
        The quality of the reconstructed tomograms strongly influences
        all downstream analyses, including particle localization,
        classification, averaging, and denoising.

        Inputs and Dataset Consistency

        The protocol requires two major inputs: aligned tilt-series data
        and corresponding CTF estimations. These datasets must be
        compatible and correctly associated through matching tilt-series
        identifiers. Proper consistency between tilt geometry and CTF
        metadata is essential because reconstruction accuracy depends on
        the correct interpretation of both spatial alignment and optical
        distortions.

        From a biological perspective, reconstruction quality depends
        heavily on the quality of the preprocessing stages. Poor tilt
        alignment, incomplete metadata, or inaccurate CTF estimation may
        generate blurred tomograms, reconstruction artifacts, or reduced
        interpretability of intracellular structures.

        Tomogram Dimensions and Sampling

        The protocol allows the user to define the reconstructed
        tomogram dimensions along X, Y, and Z. These dimensions should
        be large enough to contain the entire specimen throughout the
        reconstruction process, particularly when rotations or specimen
        tilts increase the effective occupied volume.

        Choosing the tomogram thickness is biologically important. A
        tomogram that is too thin may truncate relevant cellular
        structures or particles, while an excessively thick tomogram
        increases storage requirements and computational cost without
        necessarily improving biological interpretation. In many cases,
        selecting a thickness that approximately matches the specimen
        thickness produces the most efficient and interpretable results.

        The protocol also allows reconstruction at a larger effective
        pixel size through binning. Increasing the binned pixel size
        reduces computational demands and disk usage, which is often
        advantageous during exploratory analysis or large-scale
        processing campaigns. Smaller pixel sizes preserve more high-
        resolution information but substantially increase processing
        time and storage consumption.

        Reconstruction Modes

        The protocol supports reconstruction of either a single
        tomogram or all tomograms within the dataset. Reconstructing a
        single tomogram is useful during parameter optimization, quality
        control, or troubleshooting because it allows rapid inspection
        of reconstruction settings before processing an entire
        collection.

        Full dataset reconstruction is typically performed once the user
        has validated that alignment, sampling, and reconstruction
        parameters produce biologically meaningful results. In large
        tomography studies, beginning with a representative tomogram is
        generally recommended before scaling to the complete dataset.

        CTF Correction

        The protocol optionally performs CTF correction during the
        reconstruction process. This step compensates for optical phase
        distortions introduced during image acquisition and can improve
        interpretability of structural features within the tomogram.

        From a biological standpoint, CTF correction is particularly
        important when the reconstructed tomograms will be used for
        subtomogram averaging, template matching, or high-resolution
        refinement. However, the effectiveness of CTF correction depends
        strongly on the accuracy of the estimated defocus values and the
        quality of the tilt-series acquisition.

        In some exploratory or low-resolution workflows, users may
        initially reconstruct without CTF correction for faster
        processing and later generate corrected tomograms for detailed
        analysis.

        Odd and Even Tomograms

        The protocol can generate odd and even tomograms derived from
        separated subsets of the tilt-series data. These paired
        tomograms are particularly useful for denoising approaches such
        as cryoCARE and related machine-learning workflows.

        Biologically, odd/even tomograms preserve statistically
        independent signal representations while maintaining the same
        specimen content. This independence is essential for training
        denoising models without introducing artificial correlations or
        overfitting.

        To use this functionality effectively, the input tilt-series
        must already contain odd and even frame separation generated in
        earlier preprocessing stages.

        Tilt Geometry and Angular Offsets

        The protocol allows the introduction of a global tilt-angle
        offset during reconstruction. This option is useful in specific
        experimental geometries, such as lamella preparation workflows
        where specimens are milled at consistent angles relative to the
        microscope coordinate system.

        Correct handling of tilt geometry is biologically important
        because systematic angular offsets may otherwise distort the
        reconstructed cellular architecture or introduce anisotropic
        reconstruction artifacts.

        Outputs and Interpretation

        After completion, the protocol produces reconstructed tomograms
        ready for visualization and downstream analysis. These tomograms
        retain acquisition metadata and sampling information necessary
        for subsequent cryo-electron tomography workflows.

        When CTF correction is enabled, the resulting tomograms are
        marked accordingly and are generally more appropriate for
        subtomogram averaging and structural refinement. If odd/even
        tomograms are generated, these outputs can be directly employed
        in denoising or validation workflows.

        Biological interpretation of the reconstructed tomograms should
        always include visual inspection. Users should verify membrane
        continuity, particle visibility, reconstruction completeness,
        and the absence of strong alignment or missing-wedge artifacts.

        Practical Recommendations

        For exploratory analysis, it is often advantageous to begin with
        a moderately binned reconstruction because it significantly
        reduces processing time while preserving sufficient structural
        detail for quality assessment.

        During parameter optimization, reconstructing a single
        representative tomogram allows efficient testing of dimensions,
        tilt-angle offsets, and CTF correction strategies before
        launching large-scale processing.

        When preparing tomograms for subtomogram averaging or machine-
        learning denoising, careful selection of sampling rate and
        tomogram thickness becomes especially important. Excessively
        large tomograms increase storage and computational cost without
        necessarily improving downstream biological interpretation.

        Final Perspective

        Tomogram reconstruction is one of the foundational stages of
        cryo-electron tomography because it transforms aligned projection
        images into biologically interpretable three-dimensional
        representations of cellular organization and molecular
        architecture. Reliable reconstruction depends not only on
        computational parameters, but also on accurate preprocessing,
        appropriate sampling choices, and careful consideration of the
        biological specimen under study.
    """

    _label = 'Tomograms reconstruction'
    _possibleOutputs = tomoRecOtputObjects

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.tsDict = None
        self.ctfDict = None

        # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        form.addSection(label=Message.LABEL_INPUT)
        form.addParam('inTsSet', PointerParam,
                      pointerClass='SetOfTiltSeries',
                      important=True,
                      label='Tilt-series')
        form.addParam('inCtfSet', PointerParam,
                      pointerClass='SetOfCTFTomoSeries',
                      important=True,
                      label='CTFs')
        self._insertBinThreadsParam(form)
        form.addParam('unbinnedWidth', IntParam,
                      default=-1,
                      allowsNull=False,
                      important=True,
                      label='Unbinned tomogram width (px)',
                      help="The tomogram X-dimension in unbinned pixels. It is recommended to use "
                           "a slightly larger tomogram volume than the actual size of the images so that if "
                           "they rotate, all pixels will still be in the tomogram.\n"
                           "If set to -1, the corresponding image dimension with an extra 10% will be used, e. g., "
                           "and image with X-dim equal to 4000 px would be passed to Relion as width = 4400 px.")
        form.addParam('unbinnedHeight', IntParam,
                      default=-1,
                      allowsNull=False,
                      important=True,
                      label='Unbinned tomogram height (px)',
                      help="The tomogram Y-dimension in unbinned pixels. See the help of parameter 'Unbinned tomogram "
                           "width (px)' for more details.")
        form.addParam('unbinnedThickness', IntParam,
                      default=-1,
                      allowsNull=False,
                      important=True,
                      label='Unbinned tomogram thickness (px)',
                      help="The tomogram Z-dimension in unbinned pixels. For your own data, you may want to test a "
                           "few values here to make sure the tomogram thickness is not too small to contain your "
                           "entire sample. If you intend to denoise your tomograms later with cryoCARE, it is better "
                           "not to pick a tomogram thickness that is much greater than the thickness of your sample, "
                           "because the denoising protocol randomly extracts subtomograms from your tomograms and you "
                           "don’t want too many without signal.\n"
                           "If set to -1, a thickness of 5/16 of the image width will be passed to Relion.")
        form.addParam('binnedPixSize', FloatParam,
                      allowsNull=False,
                      important=True,
                      label='Binned pixel size (Å/px)',
                      help='The tomogram will be downscaled to this pixel size. Typically, the larger the pixel size, '
                           'the faster the tomogram reconstruction and the less space the tomograms occupy on disk.')
        form.addParam('recTomoMode', EnumParam,
                      display=EnumParam.DISPLAY_HLIST,
                      choices=['Single tomogram', 'All tomograms'],
                      default=ALL_TOMOS,
                      label='Reconstruction mode',
                      help='Choose a reconstruction option. If the option Single tomograms is selected, then the '
                            'program will only reconstruct the tomogram chosen by tilt-series id.')
        form.addParam('tomoId', StringParam,
                      condition='recTomoMode == %s' % SINGLE_TOMO,
                      label='Tomogram to be reconstructed')
        form.addParam('doCtfCorrection', BooleanParam,
                      default=False,
                      label='Perform CTF correction?')
        form.addParam('genEvenOddTomos', BooleanParam,
                      default=False,
                      label='Generate the odd/even tomograms?',
                      help="Generate the odd/even tomograms that can be used for example to denoise with cryoCARE. For "
                           "this option to work, the introduced tilt-series must have their even/odd tilt-series "
                           "generated in the motion correction step.")
        form.addParam('tiltAngleOffset', FloatParam,
                      default=0,
                      label='Tilt angle offset (deg',
                      help="The tomogram tilt angles will all be changed by this value. This may be useful to "
                           "reconstruct lamellae that are all milled under a given angle. All tomograms will be "
                           "reconstructed with the same offset.")
        # TODO: add the params related to the 2D sums of the central Z-slices?
        form.addParallelSection(threads=0, mpi=1)

    # -------------------------- INSERT steps functions -----------------------
    def _insertAllSteps(self):
        self._initialize()
        self._insertFunctionStep(self.convertInputStep, needsGPU=False)
        self._insertFunctionStep(self.reconstructTomogramsStep, needsGPU=False)
        self._insertFunctionStep(self.createOutputStep, needsGPU=False)

    # -------------------------- STEPS functions ------------------------------
    def _initialize(self):
        tsSet = self.inTsSet.get()
        ctfSet = self.inCtfSet.get()
        # Compute matching TS id among the tilt-series and the CTFs, they all could be a subset
        tsIds = tsSet.getTSIds()
        self.info("TsIds present in Tilt series are: %s" % tsIds)
        ctfTsIds = ctfSet.getTSIds()
        self.info("TsIds present in CTFs are: %s" % tsIds)
        presentTsIds = set(tsIds) & set(ctfTsIds)

        # Validate the intersection
        if len(presentTsIds) > 0:
            self.info("tsId matches between the introduces tilt-series and CTFs: %s" % presentTsIds)
        else:
            raise Exception("There isn't any common tilt-series ids among the CTFs and tilt-series introduced.")

        self.tsDict = {ts.getTsId(): ts.clone(ignoreAttrs=[]) for ts in tsSet if ts.getTsId() in presentTsIds}
        self.ctfDict = {ctf.getTsId(): ctf.clone(ignoreAttrs=[]) for ctf in ctfSet if ctf.getTsId() in presentTsIds}

    def convertInputStep(self):
        outPath = self._getExtraPath()
        writer = convert50_tomo.Writer()
        # Aligned tilt-series star files: the one corresponding to the set and each TS star file
        writer.alignedTsSet2Star(self.tsDict, outPath)
        writer.tsSet2Star(self.tsDict, self.ctfDict, outPath)

    def reconstructTomogramsStep(self):
        program = getProgram('relion_tomo_reconstruct_tomogram', nMpi=self.numberOfMpi.get())
        Plugin.runRelionTomo(self, program, self.genTomoRecCmd())

    def createOutputStep(self):
        inTsSet = self.inTsSet.get()
        outTomoSet = SetOfTomograms.create(self._getPath(), template='tomograms%s.sqlite')
        outTomoSet.copyInfo(inTsSet)
        outTomoSet.setSamplingRate(self.binnedPixSize.get())
        if self.recTomoMode.get() == SINGLE_TOMO:
            tsId = self.tomoId.get()
            ts = inTsSet.getItem(TiltSeries.TS_ID_FIELD, tsId)
            self._createTomo(ts, outTomoSet)
        else:
            for tsId, ts in self.tsDict.items():
                self._createTomo(ts, outTomoSet)

        self._defineOutputs(**{self._possibleOutputs.tomograms.name: outTomoSet})
        self._defineSourceRelation(inTsSet, outTomoSet)

    # --------------------------- UTILS functions -----------------------------
    @classmethod
    def isDisabled(cls):
        """ Return True if this Protocol is disabled.
        Disabled protocols will not be offered in the available protocols."""
        return False if IS_RELION_50 else True

    def _createTomo(self, ts, outTomoSet):
        tomo = Tomogram()
        tomo.copyInfo(ts)
        tomo.setSamplingRate(self.binnedPixSize.get())
        tomo.setFileName(self.getOutTomoFileName(ts.getTsId()))
        tomo.setOrigin()
        if self.doCtfCorrection.get():
            tomo.setCtfCorrected(True)
        outTomoSet.append(tomo)
        outTomoSet.update(tomo)

    def getOutTomoFileName(self, tsId):
        return self._getExtraPath(TOMOGRAMS_DIR, f'rec_{tsId}.mrc')

    def genTomoRecCmd(self):
        cmd = [
            f'--t {self._getExtraPath(IN_TS_STAR)}',
            f'--o {self._getExtraPath()}/',
            f'--w {self.unbinnedWidth.get()}',
            f'--h {self.unbinnedHeight.get()}',
            f'--d {self.unbinnedThickness.get()}',
            f'--binned_angpix {self.binnedPixSize.get():.3f}',
            f'--j {self.binThreads.get()}'
        ]
        if self.genEvenOddTomos.get():
            cmd.append('----generate_split_tomograms')
        if self.recTomoMode.get() == SINGLE_TOMO:
            cmd.append(f'--tn {self.tomoId.get()}')
        if self.doCtfCorrection.get():
            cmd.append('--ctf')
        return ' '.join(cmd)

    # -------------------------- INFO functions -------------------------------
    def _validate(self):
        errorMsg = []
        if self.genEvenOddTomos.get() and not self.inTsSet.get().hasOddEven():
            errorMsg.append('Odd/even tomograms cannot be generated as the introduced tilt-series do not have odd/even '
                            'tilt-series associated.')
        return errorMsg

    def _warnings(self):
        warnMsg = []
        if not self.inTsSet.get().hasAlignment():
            warnMsg.append('The introduced tilt-series seems not to have alignment data.')
        return warnMsg

    def _summary(self):
        summary = []
        if self.isFinished():
            if self.recTomoMode.get() == SINGLE_TOMO:
                summary.append('The selected tomogram was *%s*.' % self.tomoId.get())
        return summary
