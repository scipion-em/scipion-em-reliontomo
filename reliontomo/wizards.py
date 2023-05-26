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
import tkinter as tk
from tkinter import BOTH

import matplotlib as plt

from pwem.emlib.image import ImageHandler
from pyworkflow.utils import Icon

plt.use('TkAgg')
import numpy as np

from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure

from pwem.wizards import EmWizard
from pyworkflow.gui import (ListTreeProviderString, dialog, configureWeigths,
                            Dialog, RESULT_CLOSE, RESULT_YES, MessageDialog,
                            showError)
from pyworkflow.object import String
from relion.wizards import RelionWizMtfSelector
from reliontomo.protocols import (ProtRelionTomoReconstruct,
                                  ProtRelionEditParticlesStar)
from reliontomo.protocols.protocol_post_process import ProtRelionPostProcess
from tomo.objects import SetOfTomograms
from tomo.utils import getObjFromRelation
from reliontomo.protocols.protocol_prepare_data import \
    outputObjects as prepareProtOutputs

RelionWizMtfSelector._targets.append((ProtRelionPostProcess, ['mtf']))


class RelionTomoIdsWizard(EmWizard):
    tomoIdParamName = 'tomoId'
    _targets = [(ProtRelionTomoReconstruct, [tomoIdParamName])]

    def show(self, form, *args):
        relionTomoRecProt = form.protocol
        inReParticles = getattr(relionTomoRecProt.protPrepare.get(),
                                prepareProtOutputs.relionParticles.name, None)
        tomoSet = getObjFromRelation(inReParticles, relionTomoRecProt,
                                     SetOfTomograms)
        tsIds = [String(tomo.getTsId()) for tomo in tomoSet]

        # Get a data provider from the operations to be used in the tree (dialog)
        provider = ListTreeProviderString(tsIds)

        # Show the dialog
        dlg = dialog.ListDialog(form.root, "Tomograms TsIds", provider,
                                "Select one of the tomograms")

        # Set the chosen value back to the form
        form.setVar(self.tomoIdParamName, dlg.values[0].get())


class RelionWizEditParticleDisplayWindows(Dialog):
    def __init__(self, master, tittle, subtomogram, mask3D, slices, axes,
                 **args):
        self.subtomogram = subtomogram
        self.mask3D = mask3D
        self.axesToShow = axes
        self.axesLabels = list(self.axesToShow.keys())
        self.slices = slices
        self.invertContrast = True
        self.hideMask = False

        Dialog.__init__(self, master, tittle, buttons=[('Close', RESULT_CLOSE),
                                                       ('Apply', RESULT_YES)],
                        default='Apply', icons={RESULT_CLOSE: Icon.BUTTON_CLOSE,
                                                RESULT_YES: Icon.BUTTON_SELECT},
                        **args)

    def body(self, master):
        # Creating the layout where all application elements will be placed
        parentFrame = tk.Frame(master)
        parentFrame.grid(row=0, column=0, sticky='news')
        configureWeigths(parentFrame)
        self.resizable(False, False)
        self.fillWizard(parentFrame)

    def fillWizard(self, parentFrame):
        """Method to create the wizard window"""
        mainFrame = tk.Frame(parentFrame)
        mainFrame.grid(row=0, column=0, sticky='news')
        configureWeigths(mainFrame)
        mainFrame.pack()

        # Menu bar
        menubar = tk.Menu(mainFrame)
        menu_operation = tk.Menu(menubar, tearoff=False)
        menu_operation.add_command(label="Invert contrast",
                                   accelerator="(double-click)",
                                   command=self.invert_contrast
                                   )
        menu_operation.add_command(label="Hide/show mask",
                                   accelerator="(Ctrl-h)",
                                   command=self.hide_mask
                                   )
        menubar.add_cascade(menu=menu_operation, label="Operations")
        self.configure(menu=menubar)

        # Plot Frame
        plotFrame = tk.Frame(mainFrame)
        plotFrame.grid(row=0, column=0, padx=0, pady=10, sticky='news')
        configureWeigths(plotFrame)

        # Sliders Frame block
        labelFrameText = 'Shift center '
        if self.slices == 0:  # will be selected the central slice in Y
            labelFrameText += '(selected the central slice in Y Negative (Top))'
        else:  # will be selected the central slice in Z
            labelFrameText += '(selected the central slice in Z Negative (Front))'

        labelFrame = tk.LabelFrame(mainFrame, text=labelFrameText)
        labelFrame.grid(row=1, column=0, padx=0, pady=0, sticky='news')
        configureWeigths(labelFrame)

        sliderFrame = tk.Frame(labelFrame)
        sliderFrame.grid(row=1, column=0, padx=0, pady=0, sticky='news')
        fill_label = tk.Label(sliderFrame, text='')
        fill_label.grid(row=0, column=0, padx=10, pady=0, sticky='we')

        self.figure = Figure(figsize=(4, 4))
        self.ax = self.figure.add_subplot(111)
        self.ax.get_xaxis().set_visible(False)
        self.ax.get_yaxis().set_visible(False)
        self.figure.subplots_adjust(left=0, bottom=0, right=1, top=1, wspace=0,
                                    hspace=0)

        # overlapping the images
        overlay_image = self.overlay_images(self.mask3D, self.subtomogram, 0, 0)

        # Showing the images overlapped in the plot
        self.ax.imshow(overlay_image, cmap='gray')

        # Creating the canvas
        self.canvas = FigureCanvasTkAgg(self.figure, master=plotFrame)
        self.canvas.draw()
        self.canvas.get_tk_widget().grid(row=0, column=0, columnspan=2)
        self.canvas.get_tk_widget().bind("<Button-1>",
                                         self.get_sliders_position)
        self.canvas.get_tk_widget().bind("<B1-Motion>", self.update_sliders)
        self.canvas.get_tk_widget().bind("<Double-Button-1>",
                                         self.invert_contrast)
        self.canvas.get_tk_widget().bind("<Control-h>", self.hide_mask)

        # Creating the sliders
        dimension = self.subtomogram.shape[1]
        first_slider_label = tk.Label(sliderFrame, text=self.axesLabels[0])
        first_slider_label.grid(row=0, column=1, padx=10, pady=10, sticky='e')
        widgetRange = [-int(dimension / 2), int(dimension / 2)]

        # First slider
        self.first_slider = tk.Scale(sliderFrame, from_=widgetRange[0],
                                     to=widgetRange[1],
                                     orient=tk.HORIZONTAL, length=100,
                                     command=lambda
                                         value: self.update_overlay_image(
                                         self.mask3D,
                                         self.subtomogram,
                                         value,
                                         self.second_slider.get(),
                                         self.ax,
                                         self.canvas))
        self.first_slider.grid(row=0, column=2, padx=10, pady=10, sticky='we')
        self.first_slider.set(self.axesToShow[self.axesLabels[0]])

        # Second slider
        second_slider_label = tk.Label(sliderFrame, text=self.axesLabels[1])
        second_slider_label.grid(row=0, column=4, padx=10, pady=10, sticky='e')
        self.second_slider = tk.Scale(sliderFrame, from_=widgetRange[0],
                                      to=widgetRange[1],
                                      orient=tk.HORIZONTAL, length=100,
                                      command=lambda
                                          value: self.update_overlay_image(
                                          self.mask3D,
                                          self.subtomogram,
                                          self.first_slider.get(),
                                          value,
                                          self.ax,
                                          self.canvas))
        self.second_slider.grid(row=0, column=5, padx=10, pady=10, sticky='we')
        self.second_slider.set(self.axesToShow[self.axesLabels[1]])

    def get_sliders_position(self, event):
        self.last_x = int(event.x)
        self.last_y = int(event.y)
        self.first_sliderValue = self.first_slider.get()
        self.second_sliderValue = self.second_slider.get()

    def update_sliders(self, event):
        """Method to update the sliders when drag and drop the image """
        # Ratio between canvas and slider values
        ratio = 4.4
        x = event.x
        y = event.y
        moveX = -int((self.last_x - x) / ratio)
        moveY = -int((self.last_y - y) / ratio)

        self.first_slider.set(self.first_sliderValue - moveX)
        self.second_slider.set(self.second_sliderValue - moveY)

    def invert_contrast(self, event=None):
        """ Method to invert the contrast when double-click mouse action"""
        self.invertContrast = not self.invertContrast
        self.update_overlay_image(self.mask3D, self.subtomogram,
                                  self.first_slider.get(),
                                  self.second_slider.get(),
                                  self.ax, self.canvas)

    def hide_mask(self, event=None):
        """ Method to hide the mask when <ctrl-h> are preset"""
        self.hideMask = not self.hideMask
        self.update_overlay_image(self.mask3D, self.subtomogram,
                                  self.first_slider.get(),
                                  self.second_slider.get(),
                                  self.ax, self.canvas)

    def overlay_images(self, mask, subtomogram, shift_x, shift_y):
        """
        Method to overlay the images
        """
        # Normalizing the pixels values between 0 and 1
        image1_norm = mask / np.max(mask)
        image2_norm = subtomogram / np.max(subtomogram)

        # Using np.roll function to roll the images along a given axis
        shifted_image2 = np.roll(image2_norm, int(-shift_x),
                                 axis=1)
        shifted_image2 = np.roll(shifted_image2, int(-shift_y),
                                 axis=0)
        # # Invert contrast
        if self.invertContrast:
            shifted_image2 = np.subtract(1, shifted_image2)
        # Overlapping the images
        overlay_image = image1_norm * shifted_image2 if not self.hideMask else shifted_image2
        overlay_image = np.clip(overlay_image, 0, 1)

        return overlay_image

    def update_overlay_image(self, mask, subtomogram, slider_x_value,
                             slider_y_value, ax, canvas):
        """ Method to save the slider values in order to update the protocol
        parameter and draw the overlapped images """

        self.axesToShow[self.axesLabels[0]] = slider_x_value
        self.axesToShow[self.axesLabels[1]] = slider_y_value
        overlay_image = self.overlay_images(mask, subtomogram,
                                            int(slider_x_value),
                                            int(slider_y_value))
        ax.imshow(overlay_image, cmap='gray')
        canvas.draw()


class RelionWizEditParticleDisplay(EmWizard):
    _targets = [(ProtRelionEditParticlesStar, ['shiftX', 'shiftY', 'shiftZ'])]

    def readMRC(self, file):
        """Read an image using the imageHandler """
        image = ImageHandler().read(file + ':mrc')
        data = image.getData()
        return data

    def show(self, form, *args):
        relionEditParticlesStarProt = form.protocol
        averageSubTomogram = relionEditParticlesStarProt.getAverageSubTomogram()
        inMask3D = relionEditParticlesStarProt.getMask3D()

        if averageSubTomogram is not None and inMask3D is not None:
            averageSubTomogramFn = averageSubTomogram.getFileName()
            inMask3DFn = inMask3D.getFileName()
            paramName = form.wizParamName

            slices, axes = self.selectAxesToShow(paramName, form)

            # Reading the averageSubTomogram and the mask3D
            subTomogram = self.readMRC(averageSubTomogramFn)
            mask3D = self.readMRC(inMask3DFn)

            # Taking the central slice of one axis
            dimension = subTomogram.shape[2]
            central_slice = int(dimension / 2)

            if slices == 0:  # will be selected a central slice in Y
                image_data_2d = subTomogram[:, central_slice, :]
                mask = mask3D[:, central_slice, :]
            else:  # will be selected a central slice in Z
                image_data_2d = subTomogram[central_slice, :, :]
                mask = mask3D[central_slice, :, :]

            apply = RelionWizEditParticleDisplayWindows(form.getRoot(),
                                                        "Relion edit particle display wizard",
                                                        image_data_2d, mask,
                                                        slices, axes)

            if apply.result == RESULT_YES:
                if paramName == 'shiftZ':
                    form.setVar('shiftX', apply.axesToShow[apply.axesLabels[0]])
                    form.setVar('shiftZ', apply.axesToShow[apply.axesLabels[1]])
                else:
                    form.setVar('shiftX', apply.axesToShow[apply.axesLabels[0]])
                    form.setVar('shiftY', apply.axesToShow[apply.axesLabels[1]])
        else:
            showError("Input validation error",
                      "You need to select a subtomogram and a referenece mask",
                      form.getRoot())

    def selectAxesToShow(self, axis, form):
        """ Method to select the axes to show """
        label, value = self._getInputProtocol(self._targets, form.protocol)
        if axis == 'shiftZ':
            axes = {'X': value[0],
                    'Z': value[2]}  # will be selected the slice in Y
            slices = 0
        else:  # 'shiftX or shiftY'
            axes = {'X': value[0],
                    'Y': value[1]}  # will be selected the slice in Z
            slices = 1 if axis == 'shiftY' else 2
        return slices, axes
