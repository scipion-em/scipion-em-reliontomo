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
from pwem.objects import FSC, Volume
from pyworkflow.object import String
from pyworkflow.protocol import PointerParam, BooleanParam, FloatParam, GE, LE, IntParam, FileParam
from pyworkflow.utils import makePath, Message
from reliontomo import Plugin
from reliontomo.constants import POST_PROCESS_MRC, POSTPROCESS_DIR, \
    POSTPROCESS_STAR_FSC_TABLE, \
    POSTPROCESS_STAR_FSC_COLUMNS, FSC_REF_STAR, POSTPROCESS_STAR_FIELD
from reliontomo.protocols.protocol_base_relion import ProtRelionTomoBase

NO_MTF_FILE = 0


class outputObjects(Enum):
    postProcessVolume = Volume
    outputFSC = FSC


class ProtRelionPostProcess(ProtRelionTomoBase):
    """
    Sharpens refined subtomogram averaging maps and estimates masked
    Fourier Shell Correlation curves in order to improve map interpretability
    and obtain more realistic resolution measurements after refinement.

    AI Generated:

    Post-processing (ProtRelionPostProcess) — User Manual
        Overview

        The Post-processing protocol is designed to improve the visual
        and interpretative quality of reconstructed cryo-electron tomography maps
        after subtomogram refinement. In standard refinement workflows, the final
        reconstruction often appears blurred because high-resolution information
        is attenuated during imaging and reconstruction. This protocol enhances
        structural detail through sharpening procedures while also estimating
        gold-standard FSC curves that more accurately reflect the true resolution
        of the reconstructed structure.

        In practical cryo-EM and cryo-ET workflows, this step is usually
        performed immediately after a successful refinement. The resulting map is
        typically the version used for visualization, biological interpretation,
        atomic modeling, figure preparation, and deposition. The protocol is also
        important because it provides masked FSC calculations, which generally
        yield more realistic resolution estimates than unmasked curves by reducing
        the influence of solvent noise.

        Inputs and Biological Context

        The protocol requires a refined 3D volume together with its
        corresponding half maps. These half reconstructions are essential because
        gold-standard FSC estimation relies on comparing two independently refined
        datasets. Without half maps, reliable resolution assessment cannot be
        performed.

        A solvent mask is also required. This mask defines the region
        occupied by the biological structure and separates it from the surrounding
        solvent region. In cryo-EM, masking is biologically important because
        large solvent regions contribute mostly noise and can artificially lower
        FSC values. Proper masking therefore improves both sharpening behavior and
        resolution estimation.

        In most biological applications, the mask should tightly enclose
        the stable molecular region while maintaining soft edges. Soft transitions
        between protein and solvent regions help avoid Fourier artifacts and
        prevent artificially inflated resolution estimates. Excessively sharp or
        poorly designed masks may introduce distortions or misleading FSC curves.

        B-factor Sharpening

        One of the central features of this protocol is B-factor
        sharpening. Sharpening compensates for the gradual decay of high-frequency
        signal that occurs during imaging and reconstruction. The purpose is to
        enhance fine structural features such as alpha helices, beta sheets, side
        chains, membrane boundaries, or ligand densities.

        Automatic B-factor estimation is generally recommended for most
        datasets. In this mode, the protocol estimates an overall sharpening
        factor directly from the reconstruction using established cryo-EM
        approaches based on Guinier analysis. This strategy is often effective for
        well-behaved datasets with sufficient signal extending into intermediate
        resolutions.

        In some situations, however, automatic estimation may not provide
        optimal results. Highly flexible complexes, heterogeneous samples, noisy
        datasets, or maps with limited resolution range may produce unstable or
        biologically unrealistic sharpening. In these cases, users may provide a
        custom B-factor manually.

        Negative B-factors sharpen the map, whereas excessively strong
        sharpening may amplify noise and create artificial high-resolution
        features. Biological interpretation should therefore always be performed
        with caution. A map that appears visually sharper is not necessarily more
        accurate.

        FSC Weighting and Filtering

        The protocol can apply FSC-based weighting to attenuate frequency
        ranges dominated by noise. This weighting procedure improves the balance
        between signal preservation and noise suppression, often resulting in maps
        that are easier to interpret biologically.

        In some workflows, users may instead choose to disable FSC
        weighting and apply a manual low-pass filter. This approach can be useful
        when local resolution varies substantially across the structure. For
        example, rigid domains may contain high-resolution information while
        flexible regions remain poorly resolved. In such cases, a carefully chosen
        ad-hoc filter may produce maps that are visually more interpretable.

        Biological users should remember that applying overly aggressive
        filters or sharpening parameters can lead to over-interpretation of noise.
        Features near or beyond the reported resolution limit should always be
        validated carefully against known structural biology principles and, when
        possible, independent experimental evidence.

        Detector MTF Considerations

        The protocol optionally allows inclusion of the detector Modulation
        Transfer Function. Incorporating detector MTF information improves the
        physical correction of detector-dependent signal attenuation and may yield
        more accurate sharpening behavior.

        This option is particularly relevant for high-resolution datasets
        or when detector calibration is well characterized. However, many routine
        workflows still produce biologically meaningful results without explicitly
        providing detector MTF information.

        Pixel Size Calibration

        Accurate pixel size calibration is essential for meaningful FSC
        interpretation and downstream structural analysis. Small calibration
        errors may significantly affect reported resolutions, atomic fitting, and
        model validation.

        The protocol therefore allows users to provide a calibrated pixel
        size when improved measurements become available after refinement or model
        fitting. This is especially important in high-resolution studies where
        small magnification inaccuracies can influence biological conclusions.

        Outputs and Interpretation

        The protocol produces a sharpened post-processed map together with
        masked FSC curves describing the estimated resolution of the reconstruction.
        The sharpened map is generally the preferred volume for visualization and
        structural interpretation.

        The FSC curves provide quantitative information about the spatial
        frequency content supported by the data. Biological users should interpret
        these curves carefully and avoid relying exclusively on a single reported
        resolution value. Local flexibility, compositional heterogeneity, preferred
        orientations, and masking effects may all influence FSC behavior.

        Practical Recommendations

        In routine cryo-ET practice, automatic B-factor estimation with a
        carefully prepared soft solvent mask is usually the best starting point.
        Maps should always be inspected visually after sharpening to ensure that
        enhanced features remain biologically plausible.

        If the sharpened reconstruction appears noisy, fragmented, or
        contains unrealistic density features, weaker sharpening or stronger
        filtering may be appropriate. Conversely, if important structural details
        remain difficult to observe, moderate additional sharpening may improve
        interpretability.

        For flexible complexes or highly heterogeneous specimens, users
        should be particularly cautious about over-interpreting weak densities.
        Combining post-processing analysis with local resolution estimation and
        independent structural validation is often advisable.

        Final Perspective

        Post-processing is not merely a cosmetic enhancement step but a
        critical stage in the biological interpretation of cryo-EM and cryo-ET
        reconstructions. Appropriate masking, careful sharpening, and realistic
        resolution estimation strongly influence how confidently structural
        features can be interpreted. When applied thoughtfully, this protocol
        helps transform refined density maps into biologically meaningful
        representations suitable for visualization, modeling, and publication.
    """

    _label = 'Post-processing'
    _possibleOutputs = outputObjects

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        # super()._defineCommonInputParams(form)
        form.addSection(label=Message.LABEL_INPUT)
        form.addParam('inVolume', PointerParam,
                      pointerClass='Volume',
                      label='Volume to sharpen',
                      important=True,
                      help='It will provide the two unfiltered half-reconstructions that were output upon convergence '
                           'of a 3D auto-refine run.')
        form.addParam('solventMask', PointerParam,
                      pointerClass='VolumeMask',
                      important=True,
                      label='Solvent mask',
                      help='Provide a soft mask where the protein is white (1) and the solvent is black (0). '
                           'Often, the softer the mask the higher resolution estimates you will get. A soft '
                           'edge of 5-10 pixels is often a good edge width.')
        form.addParam('calPixSize', FloatParam,
                      default=-1,
                      label='Calibrated pixel size (Å/px)',
                      help='Provide the final, calibrated pixel size in Angstroms. This value may be different from '
                           'the pixel-size used thus far, e.g. when you have recalibrated the pixel size using the fit '
                           'to a PDB model. The X-axis of the output FSC plot will use this calibrated value.')
        form.addSection(label='Sharpen')
        form.addParam('estimateBFactor', BooleanParam,
                      label='Estimate B-factor automatically?',
                      default=True,
                      help='If set to Yes, then the program will use the automated procedure described by Rosenthal '
                           'and Henderson (2003, JMB) to estimate an overall B-factor for your map, and sharpen it '
                           'accordingly. Note that your map must extend well beyond the lowest resolution included '
                           'in the procedure below, which should not be set to resolutions much lower than 10 '
                           'Angstroms. ')
        form.addParam('lowestResBFit', FloatParam,
                      default=10,
                      validators=[GE(10), LE(15)],
                      condition='estimateBFactor',
                      label='Lowest resolution for auto-B factor',
                      help='This is the lowest frequency (in Angstroms) that will be included in the linear fit of '
                           'the Guinier plot as described in Rosenthal and Henderson (2003, JMB). Dont use values '
                           'much lower or higher than 10 Angstroms. If your map does not extend beyond 10 Angstroms, '
                           'then instead of the automated procedure use your own B-factor.')
        form.addParam('useOwnBFactor', BooleanParam,
                      condition='not estimateBFactor',
                      default=False,
                      label='Use your own B-factor?',
                      help='Instead of using the automated B-factor estimation, provide your own value. Use negative '
                           'values for sharpening the map. This option is useful if your map does not extend beyond '
                           'the 10A needed for the automated procedure, or when the automated procedure does not give '
                           'a suitable value (e.g. in more disordered parts of the map).')
        form.addParam('userBFactor', IntParam,
                      default=-1000,
                      validators=[LE(0)],
                      label='User-provided B-factor',
                      condition='useOwnBFactor',
                      help='Use negative values for sharpening. Be careful: if you over-sharpen your map, you may '
                           'end up interpreting noise for signal!')
        form.addParam('skipFscWeight', BooleanParam,
                      default=False,
                      label='Skip FSC-weighting?',
                      help='If set to No, then the output map will be low-pass filtered according to the '
                           'mask-corrected, gold-standard FSC-curve. Sometimes, it is also useful to provide an '
                           'ad-hoc low-pass filter, as due to local resolution variations some parts of the map may '
                           'be better and other parts may be worse than the overall resolution as measured by the FSC. '
                           'In such cases, set this option to Yes and provide an ad-hoc filter as described below.')
        form.addParam('adHocLowPassFilter', IntParam,
                      default=5,
                      validators=[GE(1)],
                      label='Ad-hoc low-pass filter (Å)',
                      condition='skipFscWeight',
                      help='This option allows one to low-pass filter the map at a user-provided frequency (in '
                           'Angstroms). When using a resolution that is higher than the gold-standard FSC-reported '
                           'resolution, take care not to interpret noise in the map for signal.')
        form.addParam('mtf', FileParam,
                      label='MTF of the detector',
                      help='If you know the MTF of your detector, provide it here. Curves for some well-known detectors'
                           ' may be downloaded from the RELION Wiki \n'
                           '- [[https://www3.mrc-lmb.cam.ac.uk/relion/index.php/ \n'
                           'Also see there for the exact format of your detector.  If you do not know the MTF of your'
                           ' detector and do not want to measure it, then by leaving this entry empty, you include '
                           'the MTF of your detector in your overall estimated B-factor upon sharpening the map.'
                           'Although that is probably slightly less accurate, the overall quality of your map '
                           'will probably not suffer very much.')
        form.addParam('origDetectorPixSize', FloatParam,
                      default=1,
                      validators=[GE(0.1), LE(2)],
                      condition='mtf',
                      label='Original detector pixel size ((Å)/pix)',
                      help='This is the original pixel size (in Angstroms) in the raw (non-super-resolution!) '
                           'micrographs.')
        self._defineExtraParams(form)

    # -------------------------- INSERT steps functions -----------------------
    def _insertAllSteps(self):
        makePath(self._getExtraPath(POSTPROCESS_DIR))
        self._insertFunctionStep(self.relionPostProcessStep, needsGPU=False)
        self._insertFunctionStep(self.createOutputStep, needsGPU=False)

    def relionPostProcessStep(self):
        Plugin.runRelionTomo(self, 'relion_postprocess', self.genPostProcessCmd())

    def createOutputStep(self):
        inVolume = self.inVolume
        fn = self._getExtraPath(FSC_REF_STAR)
        postProccesMrc = self._genPostProcessOutputMrcFile(POST_PROCESS_MRC)
        # Extend the sharpened volume with an attribute containing the postprocess.star file
        setattr(postProccesMrc, POSTPROCESS_STAR_FIELD, String(fn))

        # Output FSC
        setOfFSC = self.genFSCs(fn, POSTPROCESS_STAR_FSC_TABLE,
                                POSTPROCESS_STAR_FSC_COLUMNS)
        setattr(setOfFSC, POSTPROCESS_STAR_FIELD, String(fn))

        self._defineOutputs(**{outputObjects.postProcessVolume.name: postProccesMrc,
                               outputObjects.outputFSC.name: setOfFSC})
        self._defineSourceRelation(inVolume, postProccesMrc)

    # -------------------------- INFO functions -------------------------------
    def _validate(self):
        errorMsg = []
        if not self.inVolume.get().getHalfMaps():
            errorMsg.append('The introduced volume needs to get the corresponding half maps.')
        return errorMsg

    # --------------------------- UTILS functions -----------------------------
    def genPostProcessCmd(self):
        calPixSize = self.calPixSize.get() if self.calPixSize.get() > 0 else self.inVolume.get().getSamplingRate()
        half1, half2 = self.inVolume.get().getHalfMaps().split(',')
        cmd = ''
        cmd += '--i %s ' % half1
        cmd += '--i2 %s ' % half2
        cmd += '--o %s ' % self._getExtraPath(POSTPROCESS_DIR, POSTPROCESS_DIR.lower())
        cmd += '--mask %s ' % self.solventMask.get().getFileName()
        cmd += '--angpix %.2f ' % calPixSize
        # Sharpening
        if self.mtf.get():
            cmd += '--mtf %s ' % self.mtf.get()
            cmd += '--mtf_angpix %.2f' % self.origDetectorPixSize.get()
        if self.estimateBFactor.get():
            cmd += '--auto_bfac --autob_lowres %.2f ' % self.lowestResBFit.get()
        if self.useOwnBFactor.get():
            cmd += '--adhoc_bfac %.2f ' % self.userBFactor.get()
        # Filtering
        if self.skipFscWeight.get():
            cmd += '--skip_fsc_weighting --low_pass %i ' % self.adHocLowPassFilter.get()
        # Extra params
        cmd += self._genExtraParamsCmd()

        return cmd
