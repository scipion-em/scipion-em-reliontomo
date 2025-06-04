import logging
import re
from glob import glob
from os.path import join, exists, basename
from typing import Tuple, Union
from pyworkflow.protocol import EnumParam, NumericRangeParam
from relion.viewers import RelionViewer
from pwem.viewers import ObjectView, ChimeraView
from reliontomo.protocols import ProtRelionRefineSubtomograms, ProtRelion3DClassifySubtomograms, \
    ProtRelionDeNovoInitialModel

logger = logging.getLogger(__name__)

# Form params
VIEW_ITER = 'viewIter'
ITER_SELECT = 'iterSelection'
CLASS_MODE = 'classesMode'
CLASSES_TO_SHOW = 'classesToShow'
DISPLAY_VOL = 'displayVol'


# Viewer constants
ITER_LAST = 0
ITER_SELECTION = 1
ITER_FINAL_VOL = 2

CLASSES_ALL = 0
CLASSES_SEL = 1

VOLUME_SLICES = 0
VOLUME_CHIMERA = 1

ITER_NOTATIONS_MSG = ("Possible notation are:\n"
                      "1,5-8,10 -> [1,5,6,7,8,10]\n"
                      "2,6,9-11 -> [2,6,9,10,11]\n"
                      "2 5, 6-8 -> [2,5,6,7,8]")


class RelionTomoVolumeViewer(RelionViewer):
    """ Visualization of Relion results. """
    _targets = [ProtRelionRefineSubtomograms,
                ProtRelion3DClassifySubtomograms,
                ProtRelionDeNovoInitialModel]

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.isRefine = type(self.protocol) == ProtRelionRefineSubtomograms
        self.isCl3d = type(self.protocol) == ProtRelion3DClassifySubtomograms
        self.isInitVol = type(self.protocol) == ProtRelionDeNovoInitialModel

    def _defineParams(self, form):
        form.addSection(label='Visualization')
        viewIterChoices = ['last', 'selection']
        if not self.isCl3d:
            viewIterChoices.append('final volume')
        form.addParam(VIEW_ITER, EnumParam,
                      choices=viewIterChoices,
                      default=ITER_LAST,
                      display=EnumParam.DISPLAY_LIST,
                      label="Iteration to visualize",
                      help=f"""
        *last*: only the last iteration will be visualized.

        *selection*: you may specify a range of iterations.
        {ITER_NOTATIONS_MSG}                    

        *final volume:* depending on if the symmetry is not C1 
        and if refinement was carried out with the parameter "Run 
        in C1 and apply symmetry later" (recommended), the final 
        volume may be the symmetrized volume of the final iteration 
        once the refinement is finished."   
        """)

        form.addParam(ITER_SELECT, NumericRangeParam,
                      condition=f'viewIter == {ITER_SELECTION}',
                      label="Iterations",
                      help=f"Write the iteration list to visualize.\n"
                           f"{ITER_NOTATIONS_MSG}\n")

        group = form.addGroup('Volumes')
        multiClassCond = f'{self._getNClasses()} > 1'
        group.addParam(CLASS_MODE, EnumParam,
                       choices=['all', 'selection'],
                       default=CLASSES_ALL,
                       display=EnumParam.DISPLAY_HLIST,
                       condition=multiClassCond,
                       label='Classes to visualize')

        group.addParam(CLASSES_TO_SHOW, NumericRangeParam,
                       condition=f'{multiClassCond} & {CLASS_MODE} == {CLASSES_SEL}',
                       label='Classes',
                       help=f"Write the class list to visualize. "
                            f"If empty, all the classes will be displayed.\n"
                            f"{ITER_NOTATIONS_MSG}\n")

        group.addParam(DISPLAY_VOL, EnumParam,
                       choices=['slices', 'chimera'],
                       default=VOLUME_SLICES,
                       display=EnumParam.DISPLAY_HLIST,
                       label='Display volume with',
                       help='*slices*: display volumes as 2D slices along z axis.\n'
                            '*chimera*: display volumes as surface with Chimera.')

    def _getVisualizeDict(self) -> dict:
        visualizeDict = {DISPLAY_VOL: self._showVolumes}
        return visualizeDict

    def _getAttrib(self, attribName: str, retScipionType: bool = False):
        scipionTypeAttrib = getattr(self, attribName)
        return scipionTypeAttrib if retScipionType else scipionTypeAttrib.get()

    def _hasClasses(self):
        return self.isCl3d or self.isInitVol

    def _getInParticles(self):  # -> RelionSetOfPseudoSubtomograms:
        return self.protocol.inReParticles.get()

    def _getNClasses(self) -> int:
        return self.protocol.numberOfClasses.get() if self._hasClasses() else 1

    def _getFirstAndLastIter(self) -> Tuple[int, int]:
        filesPath = self.protocol._getExtraPath()
        filesPattern = '_half1_class001' if self.isRefine else '_class001'
        return self._readResultFiles(filesPath, filesPattern)
        # return None

    @staticmethod
    def _readResultFiles(filesPath: str,
                         filesPattern: str,
                         fileExt: str = 'mrc') -> Tuple[int, int]:
        iterations = []
        for file in glob(join(filesPath, f'*{filesPattern}*.{fileExt}')):
            # Check if the file matches the given pattern
            if re.match(fr'.*{filesPattern}\.{fileExt}$', file):
                # Extract the iteration number using regular expression
                try:
                    iteration = re.search(fr'_it(\d+){filesPattern}\.{fileExt}$', file).group(1)
                    # Append the iteration number to the list
                    iterations.append(int(iteration))
                except Exception as e:
                    logger.info(f'_readResultFiles failed for:\n'
                                f'\t- file = {file}\n'
                                f'\t- error = {e}')
                    continue

        return min(iterations), max(iterations)

    def _createVolumesSqlite(self):
        """Write a sqlite with all volumes selected for visualization."""
        path = self.protocol._getExtraPath('relion_viewer_volumes.sqlite')
        samplingRate = self._getInParticles().getSamplingRate()

        files = []
        volumes = self._getVolumeNames()
        for volFn in volumes:
            if not exists(volFn):
                raise FileNotFoundError(f"Missing volume file: {volFn}\n Please select "
                                        "a valid class or iteration number.")
            logger.debug(f"Adding vol: {volFn}")
            if not volFn.endswith(":mrc"):
                files.append(volFn + ":mrc")

        self.createVolumesSqlite(files, path, samplingRate)
        return [ObjectView(self._project, self.protocol.strId(), path)]

    def _showVolumesChimera(self):
        """ Create a chimera script to visualize selected volumes. """
        volumes = self._getVolumeNames()
        cmdFile = self.protocol._getExtraPath('chimera_volumes.cxc')
        with open(cmdFile, 'w+') as f:
            for vol in volumes:
                # We assume that the chimera script will be generated
                # at the same folder as relion volumes
                localVol = basename(vol)
                vol = vol.replace(":mrc", "")
                if exists(vol):
                    f.write("open %s\n" % localVol)
            f.write('tile\n')
        view = ChimeraView(cmdFile)
        return [view]

    def _getVolumeNames(self):
        vols = []
        nClasses = self._getNClasses()
        classesToShow = self._getAttrib(CLASSES_TO_SHOW, retScipionType=True)
        if classesToShow.get():
            selectedClasses = self._getRange(classesToShow, 'classes')
        else:
            selectedClasses = list(range(1, nClasses + 1))
        viewerIter = self._getAttrib(VIEW_ITER)
        if viewerIter == ITER_FINAL_VOL:
            logger.debug(f"Final volumes requested.")
            if self.isInitVol:
                if nClasses == 1:
                    volFn = self._getInitialModelOutFn()
                    self.__appendVolFn(volFn, vols)
                else:
                    for clId in selectedClasses:
                        volFn = self._getInitialModelOutFn(clId)
                        self.__appendVolFn(volFn, vols)
            else:  # Refinement
                volFn = self._getRefineResultFn()
                self.__appendVolFn(volFn, vols)
        else:
            # Get the iterations
            firstIter, lastIter = self._getFirstAndLastIter()
            if viewerIter == ITER_LAST:
                iterations = [lastIter]
            else:  # viewerIter == ITER_SELECTION
                iterations = self._getRange(self._getAttrib(ITER_SELECT, retScipionType=True), 'iterations')

            halfStr = '_half1' if self.isRefine else ''
            for it in iterations:
                for clId in selectedClasses:
                    volFn =self._getIterVolname(it, clId, halfStr=halfStr)
                    self.__appendVolFn(volFn, vols)

        return vols

    def _getInitialModelOutFn(self, classIndex: Union[int, None] = None) -> str:
        initialModelFn = self.protocol._getInitialModelOutFn(classIndex)
        if exists(initialModelFn):
            return initialModelFn
        else:
            raise FileNotFoundError(f"Volume {initialModelFn} does not exists.\n"
                                    "Please select a valid class number.")

    def _getIterVolname(self, iter: int, classInd: int, halfStr: str = '') -> str:
        volFn = self.protocol._getExtraPath(f'_it{iter:03d}{halfStr}_class{classInd:03d}.mrc')
        if exists(volFn):
            return volFn
        else:
            raise FileNotFoundError(f"Volume {volFn} does not exists.\n"
                                    "Please select a valid iteration and/or class number.")

    def _getRefineResultFn(self):
        volFn = self.protocol._getRefineResultFn()
        if exists(volFn):
            return volFn
        else:
            raise FileNotFoundError(f"Volume {volFn} does not exists.\n"
                                    "Please wait until the refinement finishes.")

    @staticmethod
    def __appendVolFn(volFn: str, vols: list) -> None:
        if exists(volFn):
            volFn = volFn.replace(":mrc", "")
            vols.append(volFn)
        else:
            raise FileNotFoundError(f"Volume {volFn} does not exists.\n"
                                    "Please select a valid iteration number or wait the "
                                    "protocol to finish in case of the symmetrized final volumes.")

    # -------------------------- INFO functions -------------------------------
    def _validate(self):
        pass

    # classIndexStr = f'{classIndex:03d}' if classIndex is not None else ''

# b

