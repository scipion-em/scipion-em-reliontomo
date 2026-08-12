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
from os import rename
from os.path import exists, abspath
import mrcfile
import numpy as np
from emtable import Table
from pwem.emlib.image.image_readers import Dm4ImageReader
from pyworkflow.protocol import PointerParam, IntParam, GE, BooleanParam, LEVEL_ADVANCED, FloatParam, EnumParam, \
    FileParam
from pyworkflow.utils import Message, makePath, removeBaseExt, createLink, replaceBaseExt, cyanStr
from reliontomo import Plugin
from reliontomo.constants import (IN_TS_STAR, FRAMES_DIR, MOTIONCORR_DIR,
                                  RLN_TOMO_NOMINAL_STAGE_TILT_ANGLE, RLN_MICROGRAPH_NAME, RLN_MICROGRAPH_NAME_EVEN,
                                  RLN_MICROGRAPH_NAME_ODD)
from reliontomo.convert import convert50_tomo, readTsStarFile
from reliontomo.protocols.protocol_base_relion import ProtRelionTomoBase, IS_RELION_50
from reliontomo.utils import getProgram
from tomo.objects import SetOfTiltSeries, TiltSeries

logger = logging.getLogger(__name__)

# Gain rotation values
NO_ROTATION = 0
DEG_90 = 1
DEG_180 = 2
DEF_270 = 3

# Gain flip values
NO_FLIP = 0
FLIP_UPSIDE_DOWN = 1
FLIP_LEFT_RIGHT = 2

# Suffixes
EVEN = 'even'
ODD = 'odd'

MRC = '.mrc'
MRCS = '.MRCS'


class outputObjects(Enum):
    tiltSeries = SetOfTiltSeries()
    tiltSeriesEven = SetOfTiltSeries()
    tiltSeriesOdd = SetOfTiltSeries()


class ProtRelionTomoMotionCorr(ProtRelionTomoBase):
    """
    Performs motion correction of tilt-series movies in cryo-electron
    tomography workflows in order to compensate for beam-induced specimen
    movement and generate aligned tilt images suitable for downstream
    reconstruction and subtomogram analysis.

    AI Generated:

    Motion Correction of Tilt-series Movies (ProtRelionTomoMotionCorr) — User Manual
        Overview

        The Motion Correction protocol is responsible for correcting
        beam-induced motion present in tilt-series movie data acquired during
        cryo-electron tomography experiments. During electron exposure, the
        sample often undergoes small physical displacements caused by radiation
        damage, charging effects, mechanical instability, or ice relaxation.
        These motions blur the recorded signal and reduce the attainable
        resolution of downstream reconstructions.

        The primary objective of this protocol is to align the individual
        frames of each movie so that specimen motion is compensated before further
        processing. The resulting aligned tilt images provide improved signal
        quality, sharper Thon rings for CTF estimation, and more reliable inputs
        for tomogram reconstruction and subtomogram averaging workflows.

        In biological cryo-ET studies, proper motion correction is one of
        the earliest and most critical preprocessing stages. Errors introduced at
        this point propagate through the entire workflow and may negatively affect
        alignment accuracy, contrast transfer estimation, particle picking, and
        high-resolution refinement.

        Inputs and General Workflow

        The protocol requires a collection of tilt-series movies acquired
        during tomography data collection. Each tilt image is composed of multiple
        movie frames that capture the temporal evolution of the specimen during
        electron exposure. The protocol aligns these frames to reduce motion blur
        and generates corrected tilt images organized into aligned tilt-series.

        The corrected outputs can subsequently be used for tomogram
        reconstruction, particle extraction, denoising, subtomogram refinement,
        and structural interpretation. In many cryo-ET pipelines, motion
        correction is performed immediately after data import and before any CTF
        analysis or tilt-series alignment procedures.

        Binning and Image Scaling

        The protocol allows optional binning during motion correction.
        Binning reduces the image dimensions by combining neighboring pixels and
        therefore decreases storage requirements and computational cost.

        In practical biological workflows, moderate binning is often used
        during exploratory processing or for very large super-resolution datasets.
        However, aggressive binning permanently removes high-frequency information
        and may limit the maximum achievable resolution. Users aiming for
        high-resolution subtomogram averaging generally prefer minimal binning
        during early preprocessing.

        Because binning changes the effective pixel size, downstream
        processing steps must remain consistent with the selected scaling factor.

        Patch-based Motion Correction

        Motion correction may be performed globally or locally using image
        patches. Local patch-based correction is biologically important because
        different regions of the specimen can move differently during exposure.
        This behavior is especially common in thicker cryo-ET samples, lamellae,
        cellular environments, or flexible ice regions.

        Increasing the number of patches allows the correction to model
        more spatially heterogeneous motion. However, excessive subdivision may
        reduce robustness when image contrast is weak. Smaller patch sizes are
        generally more useful for large fields of view or highly non-uniform
        motion, whereas simpler global correction may be sufficient for stable,
        thin specimens.

        Biological users should balance local flexibility against noise
        sensitivity when selecting patch parameters.

        Frame Grouping and Dose Considerations

        The protocol can group consecutive frames before estimating motion.
        Frame grouping improves signal-to-noise ratio and stabilizes alignment,
        particularly in low-dose conditions typical of tomography experiments.

        Excessive grouping, however, may reduce temporal resolution and
        limit the ability to correct rapid beam-induced movements. In practice,
        moderate grouping values often provide a good compromise between alignment
        stability and motion accuracy.

        Dose-dependent considerations are especially important in cryo-ET
        because high tilt angles typically contain lower signal and increased
        effective specimen thickness. Reliable motion correction therefore helps
        preserve weak high-resolution information across the entire tilt range.

        Gain Reference and Detector Corrections

        The protocol supports detector gain correction and detector defect
        handling. Gain references compensate for non-uniform detector response,
        while defect correction addresses permanently damaged or unreliable pixels.

        Proper detector calibration is essential for accurate intensity
        normalization and optimal image quality. Incorrect gain orientation or
        flipping may introduce artifacts that propagate throughout the tomography
        workflow.

        In biological practice, users should carefully verify detector
        configuration parameters during initial data processing, especially when
        working with new microscope setups or unfamiliar acquisition software.

        EER Data Processing

        The protocol supports Electron Event Representation data generated
        by modern direct electron detectors. EER datasets preserve very fine
        temporal sampling and provide increased flexibility during frame grouping
        and dose management.

        Fractionation settings determine how detector events are grouped
        into processing frames. Smaller fractions preserve finer temporal
        information but increase computational cost and storage demands. Larger
        fractions simplify processing and may improve robustness in lower-quality
        datasets.

        Biological users performing high-resolution tomography often
        benefit from carefully optimized EER fractionation strategies tailored to
        their specimen type and acquisition conditions.

        Even and Odd Frame Outputs

        The protocol optionally produces separate outputs from even and odd
        movie frames. These paired datasets are particularly valuable for
        downstream denoising approaches such as cryoCARE and related machine
        learning strategies.

        Although generating even and odd outputs increases processing time
        and storage requirements, many laboratories routinely enable this option
        because it preserves flexibility for future denoising and enhancement
        steps.

        Power Spectra and CTF Estimation

        Motion correction can also generate summed power spectra optimized
        for contrast transfer function estimation. Accurate CTF determination is
        especially important in tomography because signal quality decreases
        substantially at high tilt angles.

        Properly averaged power spectra improve visibility of Thon rings
        and facilitate more reliable defocus estimation. This directly influences
        subsequent tomogram reconstruction quality and subtomogram refinement
        performance.

        Outputs and Their Interpretation

        The protocol produces aligned tilt-series suitable for downstream
        tomographic workflows. Depending on the selected options, additional
        outputs may include even-frame and odd-frame tilt-series for denoising
    applications.

        The corrected images should exhibit improved sharpness, reduced
        frame blurring, and more consistent high-frequency signal. Biological
        users commonly inspect power spectra, image contrast, and visual alignment
        quality to confirm that motion correction has performed successfully.

        Practical Recommendations

        For routine cryo-ET processing, moderate patch-based alignment and
        careful gain correction generally provide robust results. Users should
        avoid excessive binning when high-resolution subtomogram averaging is the
        final objective.

        When processing challenging specimens such as thick lamellae,
        cellular tomograms, or highly tilted datasets, local motion correction and
        appropriate frame grouping become particularly important. In these cases,
        testing several parameter combinations may substantially improve final map
        quality.

        Laboratories interested in denoising or machine learning
        enhancement should strongly consider generating even and odd frame
        outputs during initial preprocessing, since recreating them later may be
        impossible without repeating the entire workflow.

        Final Perspective

        Motion correction is one of the foundational preprocessing stages
        in cryo-electron tomography. Accurate compensation of beam-induced motion
        preserves structural information that would otherwise be lost during image
        acquisition. Careful optimization of motion correction parameters can
        significantly improve tomogram quality, subtomogram averaging resolution,
        and the biological interpretability of reconstructed structures.
    """
    _label = 'Motion correction of tilt-series movies'
    _possibleOutputs = outputObjects

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.gainFn = None
        self.defectsFn = None

    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        form.addSection(label=Message.LABEL_INPUT)
        form.addParam('inputTiltSeriesM', PointerParam,
                      pointerClass='SetOfTiltSeriesM',
                      important=True,
                      label='Tilt-series movies')
        self._insertBinThreadsParam(form)
        form.addParam('outputInFloat16', BooleanParam,
                      label='Write output in float16?',
                      expertLevel=LEVEL_ADVANCED,
                      default=False,
                      help="If set to Yes, RelionCor2 will write output images in float16 MRC format. This will save a "
                           "factor of two in disk space compared to the default of writing in float32. Note that "
                           "RELION and CCPEM will read float16 images, but other programs may not (yet) do so. "
                           "Also note that this option does not work with UCSF MotionCor2.")
        form.addParam('saveEvenOdd', BooleanParam,
                      label='Save even/odd images?',
                      default=False,
                      help="It will write output images summed from both the even frames of the input movie and the "
                           "odd frames of the input movie. This generates two versions of the same movie which are "
                           "essential if you wish to carry out denoising later with cryoCARE. If you are unsure "
                           "whether you will need denoising later, it is best to select Yes, but be aware this option "
                           "increases the processing time.")
        form.addParam('saveSumPowerSpectra', BooleanParam,
                      label='Save sum of power spectra?',
                      default=True,
                      help="Sum of non-dose weighted power spectra provides better signal for CTF estimation. You "
                           "must use this option when writing in float16.")
        form.addParam('sumPowerSpectraNFrames', IntParam,
                      label='Sum power spectra every n frames',
                      default=4,
                      condition='saveSumPowerSpectra',
                      validators=[GE(0)],
                      help="McMullan et al (Ultramicroscopy, 2015) suggest summing power spectra every 4.0 e/A2 gives "
                           "optimal Thon rings.")

        form.addSection(label='Motion')
        form.addParam('bFactor', IntParam,
                      label='Bfactor',
                      default=150,
                      validators=[GE(0)],
                      help="The B-factor that will be applied to the micrographs.")
        line = form.addLine('Number of patches',
                            help='Number of patches to be used for patch based alignment.')
        line.addParam('patchX', IntParam, default=1, label='X')
        line.addParam('patchY', IntParam, default=1, label='Y')
        form.addParam('groupFrames', IntParam,
                      label='Group frames',
                      default=1,
                      validators=[GE(1)],
                      help="Average together this many frames before calculating the beam-induced shifts.")
        form.addParam('binningFactor', FloatParam,
                      label='Binning factor',
                      default=1,
                      validators=[GE(1)],
                      help="Bin the micrographs this much by a windowing operation in the Fourier Transform. Binning "
                           "at this level is hard to un-do later on, but may be useful to down-scale super-resolution "
                           "images. Float-values may be used. Do make sure though that the resulting micrograph size "
                           "is even.")
        form.addParam('gainRot', EnumParam,
                      default=NO_ROTATION,
                      choices=['No rotation (0)',
                               ' 90 degrees (1)',
                               '180 degrees (2)',
                               '270 degrees (3)'],
                      label='Gain rotation',
                      help="Rotate the gain reference by this number times 90 "
                           "degrees clockwise in relion_display. This is the "
                           "same as -RotGain in MotionCor2. \n"
                           "Note that MotionCor2 uses a different convention "
                           "for rotation so it says 'counter-clockwise'.")

        form.addParam('gainFlip', EnumParam,
                      default=NO_FLIP,
                      choices=['No flipping        (0)',
                               'Flip upside down   (1)',
                               'Flip left to right (2)'],
                      label='Gain flip',
                      help="Flip the gain reference after rotation. "
                           "This is the same as -FlipGain in MotionCor2. "
                           "0 means do nothing, 1 means flip Y (upside down) "
                           "and 2 means flip X (left to right).")

        form.addParam('defectFile', FileParam,
                      allowsNull=True,
                      label='Defects file',
                      help='Location of a UCSF MotionCor2-style defect text file or a defect map that describe the '
                           'defect pixels on the detector. Each line of a defect text file should contain four numbers '
                           'specifying x, y, width and height of a defect region. A defect map is an image (MRC or '
                           'TIFF), where 0 means good and 1 means bad pixels. The coordinate system is the same as the '
                           'input movie before application of binning, rotation and/or flipping.\n\n'
                           '_Note that the format of the defect text is DIFFERENT from the defect text produced by '
                           'SerialEM!_\n One can convert a SerialEM-style defect file into a defect map using IMOD '
                           'utilities e.g.:\n'
                           '*clip defect -D defect.txt -f tif movie.tif defect_map.tif*\n'
                           'See explanations in the SerialEM manual.\n'
                           'Leave empty if you do not have any defects, or do not want to correct for defects on your '
                           'detector.')

        form.addSection(label='EER')
        form.addParam('eerFractionation', IntParam,
                      default=8,
                      validators=[GE(1)],
                      label='EER fractionation',
                      help="The number of hardware frames to group into one fraction. This option is relevant only "
                           "for Falcon4 movies in the EER format. Note that all 'frames' in the GUI (e.g. first and "
                           "last frame for corrected sum, dose per frame) refer to fractions, not raw detector frames. "
                           "See https://www3.mrc-lmb.cam.ac.uk/relion/index.php/Image_compression#Falcon4_EER for "
                           "detailed guidance on EER processing.")

        self._defineExtraParams(form)
        form.addParallelSection(threads=0, mpi=1)

    # -------------------------- INSERT steps functions -----------------------
    def _insertAllSteps(self):
        self._insertFunctionStep(self.convertInputStep, needsGPU=False)
        self._insertFunctionStep(self.correctMotionStep, needsGPU=False)
        self._insertFunctionStep(self.createOutputStep, needsGPU=False)

    # -------------------------- STEPS functions ------------------------------
    def convertInputStep(self):
        framesPath = self._getExtraPath(FRAMES_DIR)
        makePath(framesPath)
        writer = convert50_tomo.Writer()
        writer.tsMSet2Star(self.inputTiltSeriesM.get(), self._getExtraPath())
        # Convert the gain file to MRC if it is in another format, such as dm4
        if self.gainFn:
            gainFile = self.inputTiltSeriesM.get().getGain()
            self.gainFn = self._convertCorrectionImage(gainFile)
        # The same with the defects file
        if self.defectsFn:
            defectsFn = self.defectFile.get()
            self.defectsFn = self._convertCorrectionImage(defectsFn)

    def correctMotionStep(self):
        nMpi = self.numberOfMpi.get()
        Plugin.runRelionTomo(self,
                             getProgram('relion_run_motioncorr', nMpi=nMpi),
                             self.getMotionCorrSubtomosCmd(),
                             cwd=self._getExtraPath(),
                             numberOfMpi=nMpi)

    def createOutputStep(self):
        # Create the output set
        inTsMSet = self.inputTiltSeriesM.get()
        outSRate = inTsMSet.getSamplingRate() * self.binningFactor.get()
        outTsSet = self._genOutTsSet(inTsMSet, outSRate)
        outputsDict = {outputObjects.tiltSeries.name: outTsSet}
        if self.saveEvenOdd.get():
            outTsSetEven = self._genOutTsSet(inTsMSet, outSRate, suffix=EVEN)
            outTsSetOdd = self._genOutTsSet(inTsMSet, outSRate, suffix=ODD)
            outputsDict[outputObjects.tiltSeriesEven.name] = outTsSetEven
            outputsDict[outputObjects.tiltSeriesOdd.name] = outTsSetOdd
        # Define the outputs and the relations
        self._defineOutputs(**outputsDict)
        self._defineSourceRelation(self.inputTiltSeriesM, outTsSet)
        if self.saveEvenOdd.get():
            self._defineSourceRelation(self.inputTiltSeriesM, outTsSetEven)
            self._defineSourceRelation(self.inputTiltSeriesM, outTsSetOdd)

    def _genOutTsSet(self, inTsMSet, outSRate, suffix=''):
        outTsSet = SetOfTiltSeries.create(self._getPath(), template='tiltseries', suffix=suffix)
        outTsSet.copyInfo(inTsMSet)
        outTsSet.setSamplingRate(outSRate)
        isEvenOdd = True if suffix else False
        # Fill it with the generated tilt-series
        for tsM in inTsMSet:
            tsId = tsM.getTsId()
            outTsStarName = self.getOutTsStarFileName(tsId)
            # Rename each TS output star files as they preserve the same base name as the input files, which are
            # preceded by an in_ suffix to avoid confusion. Only for the complete TS
            if not suffix and not exists(outTsStarName):
                rename(self.getOutTsStarFileName(tsId, preffix='in'), outTsStarName)
            newTs = TiltSeries(tsId=tsId)
            outTsSet.append(newTs)
            newTs.copyInfo(tsM)
            newTs.setSamplingRate(outSRate)
            if not suffix:
                self.mountStack(newTs)  # It mounts also the even/odd if requested
            readTsStarFile(tsM, newTs, outTsStarName, self.getOutStackName(tsId, suffix=suffix),
                           self._getExtraPath(), isEvenOdd=isEvenOdd)
            outTsSet.update(newTs)
        return outTsSet

    # -------------------------- INFO functions -------------------------------
    def _validate(self):
        errorMsg = []
        if self.outputInFloat16.get() and not self.saveSumPowerSpectra.get():
            errorMsg.append("'Save sum of power spectra?' must be set to Yes when writing in float16.")
        return errorMsg

    # --------------------------- UTILS functions -----------------------------
    @classmethod
    def isDisabled(cls):
        """ Return True if this Protocol is disabled.
        Disabled protocols will not be offered in the available protocols."""
        return False if IS_RELION_50 else True

    def getMotionCorrSubtomosCmd(self):
        cmd = f'--use_own --j {self.binThreads.get()} '
        cmd += f'--i {IN_TS_STAR} '
        cmd += f'--o {MOTIONCORR_DIR}/ '
        if self.outputInFloat16.get():
            cmd += '--float16 '
        if self.saveEvenOdd.get():
            cmd += '--even_odd_split '
        if self.saveSumPowerSpectra.get():
            cmd += f'--grouping_for_ps {self.sumPowerSpectraNFrames.get()} '
        # MOTION TAB
        cmd += f'--bfactor {self.bFactor.get()} '
        cmd += f'--patch_x {self.patchX.get()} '
        cmd += f'--patch_y {self.patchY.get()} '
        cmd += f'--group_frames {self.groupFrames.get()} '
        cmd += f'--bin_factor {self.binningFactor.get()} '
        if self.gainFn:
            cmd += f'--gainref {self.gainFile} '
            cmd += f'--gain_rot {self.gainRot.get()} '
            cmd += f'--gain_flip {self.gainFlip.get()} '
        if self.defectsFn:
            cmd += f'--defect_file {self.defectsFn.get()} '
        cmd += f'--eer_grouping {self.eerFractionation.get()}'
        # EXTRA PARAMS
        cmd += self._genExtraParamsCmd()
        return cmd

    def getOutTsStarFileName(self, tsId, preffix=''):
        bName = f'{preffix}_{tsId}' if preffix else tsId
        return self._getExtraPath(MOTIONCORR_DIR, bName + '.star')

    def getOutStackName(self, tsId, suffix=''):
        bName = f'{tsId}_{suffix}' if suffix else tsId
        return self._getExtraPath(MOTIONCORR_DIR, bName + MRCS)

    def mountStack(self, ts):
        tsId = ts.getTsId()
        sRate = ts.getSamplingRate()
        dataTable = Table()
        dataTable.read(self.getOutTsStarFileName(tsId), tableName=tsId)
        dataTable.sort(RLN_TOMO_NOMINAL_STAGE_TILT_ANGLE)  # Sort by tilt angle
        self._mountCurrentStack(tsId, sRate, dataTable)
        # Mount the even/odd stacks if requested
        if self.saveEvenOdd.get():
            self._mountCurrentStack(tsId, sRate, dataTable, imgField=RLN_MICROGRAPH_NAME_EVEN, suffix=EVEN)
            self._mountCurrentStack(tsId, sRate, dataTable, imgField=RLN_MICROGRAPH_NAME_ODD, suffix=ODD)

    def _mountCurrentStack(self, tsId, sRate, dataTable, imgField=RLN_MICROGRAPH_NAME, suffix=''):
        outStackFile = self.getOutStackName(tsId, suffix=suffix)
        logger.info(f'Mounting the stack file {outStackFile}')
        alignedImgs = [self._getExtraPath(row.get(imgField)) for row in dataTable]

        # Read the first image to get the dimensions
        with mrcfile.mmap(alignedImgs[0], mode='r+') as mrc:
            img = mrc.data
            nx, ny = img.shape

        # Create an empty array in which the stack of images will be stored
        shape = (len(alignedImgs), nx, ny)
        stackArray = np.empty(shape, dtype=img.dtype)

        # Fill it with the images sorted by angle
        for i, img in enumerate(alignedImgs):
            with mrcfile.mmap(img) as mrc:
                logger.info(f'Inserting image - index [{i}], {img}')
                stackArray[i] = mrc.data

        # Save the stack in a new mrc file
        with mrcfile.new_mmap(outStackFile, shape, overwrite=True) as mrc:
            mrc.set_data(stackArray)
            mrc.update_header_from_data()
            mrc.update_header_stats()
            mrc.voxel_size = sRate

    def _getGainMrcFName(self, inGainFileName):
        return self._getExtraPath(removeBaseExt(inGainFileName) + MRC)

    def _convertCorrectionImage(self, corrImgFn: str) -> str:
        """ Reimplement to convert dm4 gain only. """
        if corrImgFn.endswith(".dm4"):
            # Get final correction image file
            finalName = self._getExtraPath(replaceBaseExt(corrImgFn, "mrc"))

            if not exists(finalName):
                Dm4ImageReader.dmToMrc(corrImgFn, finalName)
                logger.info(cyanStr(f"Converting {corrImgFn} to {finalName}"))

            return abspath(finalName)
        else:
            return corrImgFn
