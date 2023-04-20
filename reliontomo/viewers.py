from pwem.viewers import DataViewer, MODE, MODE_MD, VISIBLE
from reliontomo.objects import RelionSetOfPseudoSubtomograms
DataViewer.registerConfig(RelionSetOfPseudoSubtomograms,
                          config={MODE: MODE_MD,
                                  VISIBLE: 'id _filename _volName _coordinate._x _coordinate._y _coordinate._z '
                                           '_transform._matrix '})
