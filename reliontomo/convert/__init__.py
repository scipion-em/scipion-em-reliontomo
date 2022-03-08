from emtable import Table

from reliontomo import Plugin
from reliontomo.constants import TOMO_NAME_30
from reliontomo.convert import convert40_tomo, convert30_tomo


def createWriterTomo30(**kwargs):
    """ Create a new Writer instance for Relion 3."""
    Writer = convert30_tomo.Writer
    return Writer(**kwargs)


def createWriterTomo40(**kwargs):
    """ Create a new Writer instance for Relion 4."""
    if Plugin.isReconstringParticleFromSubtomos():
        Writer = convert30_tomo.Writer
    else:
        Writer = convert40_tomo.Writer

    return Writer(**kwargs)


def writeSetOfPseudoSubtomograms(imgSet, starFile, **kwargs):
    return createWriterTomo40(**kwargs).pseudoSubtomograms2Star(imgSet, starFile)


def writeSetOfSubtomograms(coordSet, starFile, **kwargs):
    """ Convenience function to write a SetOfSubtomograms as Relion metadata using a Writer."""
    if Plugin.isRe40():  # TODO: check the file labels to avoid plugin version checking as re3 is expected to be deprecated in Scipion, but retrocomp may be kept
        return createWriterTomo40(**kwargs).coordinates2Star(coordSet, starFile, **kwargs)
    else:
        return createWriterTomo30(**kwargs).writeSetOfSubtomograms(coordSet, starFile, **kwargs)


def writeSetOfTomograms(imgSet, starFile, **kwargs):
    """ Convenience function to write a SetOfTomograms as Relion metadata using a Writer."""
    # There's no tomogram star file in Relion3
    return createWriterTomo40(**kwargs).tiltSeries2Star(imgSet, starFile, **kwargs)


# def createReaderTomo30(starFile=None, **kwargs):
#     """ Create a new Reader instance for Relion 3."""
#     reader = convert30_tomo.Reader(**kwargs)
#     if starFile:
#         reader.read(starFile)
#     return reader
#
#
# def createReaderTomo40(**kwargs):
#     """ Create a new Reader instance for Relion 4."""
#     Reader = convert40_tomo.Reader
#     return Reader(**kwargs)


def createReaderTomo(starFile=None, **kwargs):
    if starFile:
        dataTable = Table()
        dataTable.read(starFile)
        labels = dataTable.getColumnNames()
        if TOMO_NAME_30 in labels:
            reader = convert30_tomo.Reader(**kwargs)
            isReader40 = False
        else:
            reader = convert40_tomo.Reader(**kwargs)
            isReader40 = True
        reader.read(starFile)
    else:
        if Plugin.isRe40():
            reader = convert40_tomo.Reader(**kwargs)
            isReader40 = True
        else:
            reader = convert30_tomo.Reader(**kwargs)
            isReader40 = False

    return reader, isReader40

# def readSetOfCoordinates3D(starFile, prot):
#     """ Convenience function to write a SetOfSubtomograms as Relion metadata using a Reader."""
#     # SetOfPseudoSubtomograms don't exist in Relion3
#     reader = createReaderTomo30(starFile=starFile)
#     return reader.starFile2Coords3D(prot)


def readSetOfPseudoSubtomograms(starFile, outputSet, precedentSet, vTomoScaleFactor):
    """ Convenience function to write a SetOfPseudoSubtomograms as Relion metadata using a Reader."""
    # Subtomogras are represented in Relion 4 as Pseudosubtomograms
    reader, _ = createReaderTomo()
    return reader.starFile2PseudoSubtomograms(starFile, outputSet, precedentSet, vTomoScaleFactor)


# class ClassesLoader:
#     """ Helper class to read classes information from star files produced
#     by Relion classification runs (2D or 3D).
#     """
#     def __init__(self, protocol, alignType):
#         self._protocol = protocol
#         self._alignType = alignType
#
#     def _loadClassesInfo(self, iteration):
#         """ Read some information about the produced Relion 3D classes
#         from the *model.star file.
#         """
#         self._classesInfo = {}  # store classes info, indexed by class id
#
#         modelFn = self._protocol._getFileName('model', iter=iteration)
#         modelIter = Table.iterRows('model_classes@' + modelFn)
#
#         for classNumber, row in enumerate(modelIter):
#             index, fn = relionToLocation(row.rlnReferenceImage)
#             # Store info indexed by id
#             self._classesInfo[classNumber + 1] = (index, fn, row)
#
#     def fillClassesFromIter(self, clsSet, iteration):
#         """ Create the SetOfClasses3D from a given iteration. """
#         prot = self._protocol  # shortcut
#         self._loadClassesInfo(iteration)
#
#         tableName = 'particles@'  # if Plugin.IS_GT30() else ''
#         dataStar = prot._getFileName('data', iter=iteration)
#
#         pixelSize = prot.inputParticles.get().getSamplingRate()
#         self._reader = createReader(alignType=self._alignType,
#                                     pixelSize=pixelSize)
#
#         mdIter = Table.iterRows(tableName + dataStar, key='rlnImageId')
#         clsSet.classifyItems(updateItemCallback=self._updateParticle,
#                              updateClassCallback=self._updateClass,
#                              itemDataIterator=mdIter,
#                              doClone=False)
#
#     def _updateParticle(self, item, row):
#         item.setClassId(row.rlnClassNumber)
#         self._reader.setParticleTransform(item, row)
#
#         # Try to create extra objects only once if item is reused
#         if not hasattr(item, '_rlnNormCorrection'):
#             item._rlnNormCorrection = Float()
#             item._rlnLogLikeliContribution = Float()
#             item._rlnMaxValueProbDistribution = Float()
#
#         item._rlnNormCorrection.set(row.rlnNormCorrection)
#         item._rlnLogLikeliContribution.set(row.rlnLogLikeliContribution)
#         item._rlnMaxValueProbDistribution.set(row.rlnMaxValueProbDistribution)
#
#         if hasattr(item, '_rlnGroupName'):
#             item._rlnGroupName.set(row.rlnGroupName)
#         elif hasattr(row, 'rlnGroupName'):
#             item._rlnGroupName = String(row.rlnGroupName)
#
#     def _updateClass(self, item):
#         classId = item.getObjId()
#         if classId in self._classesInfo:
#             index, fn, row = self._classesInfo[classId]
#             item.setAlignment(self._alignType)
#             if self._alignType == ALIGN_PROJ:
#                 fn += ':mrc'  # mark reference as a MRC volume
#             item.getRepresentative().setLocation(index, fn)
#             item._rlnclassDistribution = Float(row.rlnClassDistribution)
#             item._rlnAccuracyRotations = Float(row.rlnAccuracyRotations)
#             item._rlnAccuracyTranslationsAngst = Float(row.rlnAccuracyTranslationsAngst)
#             # if Plugin.IS_GT30():
#             #     item._rlnAccuracyTranslationsAngst = Float(row.rlnAccuracyTranslationsAngst)
#             # else:
#             #     item._rlnAccuracyTranslations = Float(row.rlnAccuracyTranslations)
