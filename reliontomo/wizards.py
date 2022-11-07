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
from pwem.wizards import EmWizard
from pyworkflow.gui import ListTreeProviderString, dialog
from pyworkflow.object import String
from relion.wizards import RelionWizMtfSelector
from reliontomo.protocols import ProtRelionTomoReconstruct
from reliontomo.protocols.protocol_post_process import ProtRelionPostProcess
from tomo.objects import SetOfTomograms
from tomo.utils import getObjFromRelation
from reliontomo.protocols.protocol_prepare_data import outputObjects as prepareProtOutputs

RelionWizMtfSelector._targets.append((ProtRelionPostProcess, ['mtf']))


class RelionTomoIdsWizard(EmWizard):

    tomoIdParamName = 'tomoId'
    _targets = [(ProtRelionTomoReconstruct, [tomoIdParamName])]

    def show(self, form):
        relionTomoRecProt = form.protocol
        inReParticles = getattr(relionTomoRecProt.protPrepare.get(), prepareProtOutputs.relionParticles.name, None)
        tomoSet = getObjFromRelation(inReParticles, relionTomoRecProt, SetOfTomograms)
        tsIds = [String(tomo.getTsId()) for tomo in tomoSet]

        # Get a data provider from the operations to be used in the tree (dialog)
        provider = ListTreeProviderString(tsIds)

        # Show the dialog
        dlg = dialog.ListDialog(form.root, "Tomograms TsIds", provider,
                                "Select one of the tomograms")

        # Set the chosen value back to the form
        form.setVar(self.tomoIdParamName, dlg.values[0].get())


