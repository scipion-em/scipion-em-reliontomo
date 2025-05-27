import logging
import re
from glob import glob
from os.path import join, exists, basename
from typing import Tuple, Union
from emtable import Table
from pyworkflow.protocol import EnumParam, NumericRangeParam
from relion.viewers import RelionViewer, RelionPlotter, protected_show
from pwem.viewers import ObjectView, ChimeraView
# from reliontomo.constants import RLN_CLASSNUMBER, RLN_ANGLEROT, RLN_ANGLETILT, PARTICLES_TABLE
# from reliontomo.objects import RelionSetOfPseudoSubtomograms
from reliontomo.protocols import ProtRelionRefineSubtomograms, ProtRelion3DClassifySubtomograms, \
    ProtRelionDeNovoInitialModel

logger = logging.getLogger(__name__)

# Viewer constants
ITER_LAST = 0
ITER_SELECTION = 1
ITER_SYMMETRIZED = 2

CLASSES_ALL = 0
CLASSES_SEL = 1

VOLUME_SLICES = 0
VOLUME_CHIMERA = 1

ANGDIST_2DPLOT = 0
ANGDIST_CHIMERA = 1


# DataViewer.registerConfig(RelionSetOfPseudoSubtomograms,
#                           config={MODE: MODE_MD,
#                                   VISIBLE: 'id _filename _volName _coordinate._x _coordinate._y _coordinate._z '
#                                            '_transform._matrix '})


class RelionTomoViewer(RelionViewer):
    """ Visualization of Relion results. """
    _targets = [ProtRelionRefineSubtomograms, ProtRelion3DClassifySubtomograms,
                ProtRelionDeNovoInitialModel]

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.isRefine = type(self.protocol) == ProtRelionRefineSubtomograms
        self.isCl3d = type(self.protocol) == ProtRelion3DClassifySubtomograms
        self.isInitVol = type(self.protocol) == ProtRelionDeNovoInitialModel

    def _defineParams(self, form):
        form.addSection(label='Visualization')
        form.addParam('viewIter', EnumParam,
                      choices=['last', 'selection', 'final volume'],
                      default=ITER_LAST,
                      display=EnumParam.DISPLAY_LIST,
                      label="Iteration to visualize",
                      help="""
        *last*: only the last iteration will be visualized.
        
        *selection*: you may specify a range of iterations.
        Examples:
        "1,5-8,10" -> [1,5,6,7,8,10]
        "2,6,9-11" -> [2,6,9,10,11]
        "2 5, 6-8" -> [2,5,6,7,8]                      
                                   
        *final volume:* depending on if the symmetry is not C1 
        and if refinement was carried out with the parameter "Run 
        in C1 and apply symmetry later" (recommended), the final 
        volume may be the symmetrized volume of the final iteration 
        once the refinement is finished."   
        """)
        form.addParam('iterSelection', NumericRangeParam,
                      condition='viewIter==%d' % ITER_SELECTION,
                      label="Iterations list",
                      help="Write the iteration list to visualize.")

        group = form.addGroup('Volumes')

        # if self._hasClasses():
        #     group.addParam('showClasses3D', EnumParam,
        #                    default=CLASSES_ALL,
        #                    choices=['all', 'selection'],
        #                    display=EnumParam.DISPLAY_HLIST,
        #                    label='3D Class to visualize',
        #                    help='')
        #     group.addParam('class3DSelection', NumericRangeParam,
        #                    default='1',
        #                    condition='showClasses3D == %d' % CLASSES_SEL,
        #                    label='Classes list',
        #                    help='')
        # else:
        #     if self.isInitVol:
        #         group.addHidden('showHalves', IntParam, default=3)
        #     else:
        #         group.addParam('showHalves', EnumParam, default=0,
        #                        choices=['half1', 'half2', 'both', 'final'],
        #                        label='Volume to visualize',
        #                        help='Select which half do you want to visualize.')

        group.addParam('displayVol', EnumParam,
                       choices=['slices', 'chimera'],
                       default=VOLUME_SLICES,
                       display=EnumParam.DISPLAY_HLIST,
                       label='Display volume with',
                       help='*slices*: display volumes as 2D slices along z axis.\n'
                            '*chimera*: display volumes as surface with Chimera.')
        group.addParam('displayAngDist', EnumParam,
                       choices=['2D plot', 'chimera'],
                       default=ANGDIST_2DPLOT,
                       display=EnumParam.DISPLAY_HLIST,
                       label='Display angular distribution',
                       help='*2D plot*: display angular distribution as interactive 2D in matplotlib.\n'
                            '*chimera*: display angular distribution using Chimera with red spheres.')
        # group.addParam('spheresScale', IntParam,
        #                default=-1,
        #                expertLevel=LEVEL_ADVANCED,
        #                condition='displayAngDist == %d' % ANGDIST_CHIMERA,
        #                label='Spheres distance',
        #                help='If the value is -1 then the distance is set '
        #                     'to 0.75 * xVolDim')

    def _getVisualizeDict(self) -> dict:
        # self._load()
        visualizeDict = {'displayVol': self._showVolumes,
                         # 'showClassesOnly': self._showClassesOnly,
                         'displayAngDist': self._showAngularDistribution,}
        return visualizeDict

    # def _load(self):
    #     from matplotlib.ticker import FuncFormatter
    #     self._plotFormatter = FuncFormatter(self._formatFreq)

    def _hasClasses(self):
        return self.isCl3d or self.isInitVol

    def _getInParticles(self):# -> RelionSetOfPseudoSubtomograms:
        return self.protocol.inReParticles.get()

    def _getNClasses(self) -> int:
        return self.protocol.numberOfClasses.get()

    def _getFirstAndLastIter(self) -> Union[Tuple[int, int], None]:
        filesPath = self.protocol._getExtraPath()
        if self.isInitVol:
            filesPattern = '_class001'
            return self._readResultFiles(filesPath, filesPattern)
        return None

    @staticmethod
    def _readResultFiles(filesPath: str,
                         filesPattern: str,
                         fileExt: str = 'mrc') -> Tuple[int, int]:
        iterations = []
        for file in glob(join(filesPath, f'*{filesPattern}*.{fileExt}')):
            # Check if the file matches the given pattern
            if re.match(fr'.*{filesPattern}\.{fileExt}$', file):
                # Extract the iteration number using regular expression
                iteration = re.search(fr'_it(\d+){filesPattern}\.{fileExt}$', file).group(1)
                # Append the iteration number to the list
                iterations.append(int(iteration))

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
        viewerIter = self.viewIter.get()
        if viewerIter == ITER_SYMMETRIZED:
            logger.debug(f"Symmetrized final volumes requested.")
            if nClasses == 1:
                volFn = self.protocol._getInitialModelOutFn()
                self.__appendVolFn(volFn, vols)
            else:
                for clId in range(nClasses):
                    volFn = self.protocol._getInitialModelOutFn(clId + 1)
                    self.__appendVolFn(volFn, vols)
        else:
            # Get the iterations
            firstIter, lastIter = self._getFirstAndLastIter()
            if viewerIter == ITER_LAST:
                iterations = [lastIter]
            elif viewerIter == ITER_SELECTION:
                iterations = self._getRange(self.iterSelection, 'iterations')

            logger.debug(f"self._iterations: {iterations}.")
            for it in iterations:
                for clId in range(nClasses):
                    volFn = self.protocol._getExtraPath(f'_it{it:03d}_class{clId + 1:03d}.mrc')
                    self.__appendVolFn(volFn, vols)

        return vols

    @staticmethod
    def __appendVolFn(volFn: str, vols: list) -> None:
        if exists(volFn):
            volFn = volFn.replace(":mrc", "")
            vols.append(volFn)
        else:
            raise FileNotFoundError(f"Volume {volFn} does not exists.\n"
                                    "Please select a valid iteration number or wait the "
                                    "protocol to finish in case of the symmetrized final volumes.")

    # =============================================================================
    # showAngularDistribution
    # =============================================================================
    @protected_show
    def _showAngularDistribution(self, paramName=None):
        views = []
        _, lastIter = self._getFirstAndLastIter()
        if self.displayAngDist == ANGDIST_CHIMERA:
            for it in self._iterations:
                views.append(self._createAngDistChimera(it))

        elif self.displayAngDist == ANGDIST_2DPLOT:
            plot = self._createAngDist2D(lastIter)
            views.append(plot)
        return views

    def _createAngDist2D(self, it: int) -> RelionPlotter:
        # Common variables to use
        nparts = self._getInParticles().getSize()
        # prefixes = self._getPrefixes()
        # nrefs = len(self._refsList)
        # n = nrefs * len(prefixes)
        nClasses = self._getNClasses()
        gridsize = self._getGridSize(nClasses)

        if self.viewIter.get() == ITER_LAST:
            title = "Final"
        else:
            title = 'Iteration %d' % it

        plotter = RelionPlotter(x=gridsize[0], y=gridsize[1],
                                mainTitle=title, windowTitle="Angular Distribution")

        for clId in range(nClasses):
            mdList = self._readParticlesStarFile(it, clId)
            classId = clId + 1
            title = f'class {classId:03d}'
            sqliteFn = self.protocol._getExtraPath(f'class{classId:03d}_projections.sqlite')
            if not exists(sqliteFn):
                self.createAngDistributionSqlite(sqliteFn,
                                                 nparts,
                                                 itemDataIterator=self._iterAngles(mdList))
            plotter.plotAngularDistributionFromMd(sqliteFn, title)

        return plotter


    @staticmethod
    def _getStarFile(iteration: int) -> str:
        return f'_it{iteration:03d}_data.star'

    def _readParticlesStarFile(self, iteration: int, classIndex: int) -> list:
        from reliontomo.constants import RLN_CLASSNUMBER, PARTICLES_TABLE
        starFile = self._getStarFile(iteration)
        dataTable = Table(fileName=starFile, tableName=PARTICLES_TABLE)
        mdList = []
        for row in dataTable:
            if int(row.get(RLN_CLASSNUMBER, -1)) == classIndex:
                mdList.append(row)
        return mdList

    def _iterAngles(self, mdList: list):
        from reliontomo.constants import RLN_ANGLEROT, RLN_ANGLETILT
        for row in mdList:
            rot = float(row.get(RLN_ANGLEROT, 0))
            tilt = float(row.get(RLN_ANGLETILT, 0))
            yield rot, tilt

    # -------------------------- INFO functions -------------------------------
    def _validate(self):
        pass

    # classIndexStr = f'{classIndex:03d}' if classIndex is not None else ''


# b

