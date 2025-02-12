#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 11 20:04:58 2025

@author: Thiago Bezerra
"""

import sys
import uproot
import numpy as np
import awkward as ak
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from PyQt5.QtWidgets import QApplication, QMainWindow, QVBoxLayout, QWidget, QPushButton
from PyQt5.QtWidgets import QHBoxLayout, QSpinBox, QLabel, QRadioButton, QGroupBox

class EventDisplay(QMainWindow):
    def __init__(self, input_file):
        super().__init__()
        self.input_file = input_file
        self.initUI()
        self.colorbar2 = None

    def initUI(self):
        self.setWindowTitle('Event Display')
        self.setGeometry(100, 100, 1200, 600)  # Adjusted window size

        central_widget = QWidget()
        self.setCentralWidget(central_widget)
        layout = QVBoxLayout(central_widget)

        hBox = QHBoxLayout()

        # MC group box
        group_MC_box = QGroupBox("MC")
        vbox1 = QVBoxLayout()        
        self.signal_radio = QRadioButton('nPE')
        self.time_radio = QRadioButton('First Hit Time')
        self.med_time_radio = QRadioButton('Median Hit Time')
        vbox1.addWidget(self.signal_radio)
        vbox1.addWidget(self.time_radio)
        vbox1.addWidget(self.med_time_radio)
        group_MC_box.setLayout(vbox1)
        
        self.signal_radio.setChecked(True)  # Set default selection
        self.med_time_radio.setDisabled(True) # functionality not available yet

        # DATA group box
        group_Data_box = QGroupBox("Data")
        vbox2 = QVBoxLayout()        
        self.charge_radio = QRadioButton('SiPM Charge')
        self.rise_time_radio = QRadioButton('SiPM Time')
        vbox2.addWidget(self.charge_radio)
        vbox2.addWidget(self.rise_time_radio)
        group_Data_box.setLayout(vbox2)

        self.charge_radio.setDisabled(True)
        self.rise_time_radio.setDisabled(True)

        hBox.addWidget(group_MC_box)
        hBox.addWidget(group_Data_box)

        layout.addLayout(hBox)

        # Create a layout for the histograms
        self.hist_layout = QHBoxLayout()
        layout.addLayout(self.hist_layout)

        # Create figures for the histograms
        self.figure1, self.ax1 = plt.subplots()
        self.figure2, self.ax2 = plt.subplots()

        # Add the figures to the layout
        self.hist_layout.addWidget(self.figure1.canvas)
        self.hist_layout.addWidget(self.figure2.canvas)

        control_layout = QHBoxLayout()
        button_layout  = QVBoxLayout()

        # Event selector
        self.entry_label = QLabel('Select Entry:')
        button_layout.addWidget(self.entry_label)
        self.entry_spinbox = QSpinBox()
        self.entry_spinbox.setMinimum(0)
        button_layout.addWidget(self.entry_spinbox)

        # Display button
        self.plot_button = QPushButton('Display event')
        self.plot_button.clicked.connect(self.plot_data)
        button_layout.addWidget(self.plot_button)

        control_layout.addLayout(button_layout)
        
        self.figure3, self.ax3 = plt.subplots()
        control_layout.addWidget(self.figure3.canvas)
        
        layout.addLayout(control_layout)

    def plot_data(self):
        try:
            # Read data from ROOT file
            file = uproot.open(self.input_file)
            tree = file['output']
            entry_index = self.entry_spinbox.value()

            x = tree['mcPEx'].array(entry_start=entry_index, entry_stop=entry_index+1)[0]        
            x = ak.to_numpy(x)

            y = tree['mcPEy'].array(entry_start=entry_index, entry_stop=entry_index+1)[0]        
            y = ak.to_numpy(y)

            z = tree['mcPEz'].array(entry_start=entry_index, entry_stop=entry_index+1)[0]
            z = ak.to_numpy(z)

            t = tree['mcPEHitTime'].array(entry_start=entry_index, entry_stop=entry_index+1)[0]
            t = ak.to_numpy(t)

            # get back and front channels
            negative_pez_indices = z < 0
            positive_pez_indices = z >= 0

            x_back = x[negative_pez_indices]
            y_back = y[negative_pez_indices]
            t_back = t[negative_pez_indices]

            x_front = x[positive_pez_indices]
            y_front = y[positive_pez_indices]
            t_front = t[positive_pez_indices]

            # Define fixed ranges for the histograms
            x_range = (-1000, 1000)
            y_range = (-1000, 1000)

            if self.signal_radio.isChecked():
                # Signal Strength data handling
                print("Signal Strength selected")

                # Create 2D histograms
                h_back, xedges1, yedges1 = np.histogram2d(x_back, y_back, bins=100,range=[x_range, y_range])
                h_front, xedges2, yedges2 = np.histogram2d(x_front, y_front, bins=100,range=[x_range, y_range])            


            elif self.time_radio.isChecked():
                # Time data handling
                print("First Hit Time selected")
                
                bins = (100,100)
                # Define the edges for the bins
                xedges1 = np.linspace(x_range[0], x_range[1], bins[0] + 1)
                yedges1 = np.linspace(y_range[0], y_range[1], bins[1] + 1)      
                xedges2 = xedges1
                yedges2 = yedges1

                hit_dict_f = {}
                for xi, yi, ti in zip(x_front,y_front,t_front): #loop over vector length
                    if (xi, yi) not in hit_dict_f: #getting the hit times of each fibre
                        hit_dict_f[(xi, yi)] = ti
                    else:
                        hit_dict_f[(xi, yi)] = min(hit_dict_f[(xi, yi)], ti) 

                h_front = np.zeros(bins)
                # Fill the histogram with values from the dictionary
                for (xpos, ypos), tiso in hit_dict_f.items():
                    xi = np.digitize(xpos, xedges1) - 1
                    yi = np.digitize(ypos, yedges1) - 1
                    h_front[xi, yi] = tiso
                

                hit_dict_b = {}
                for xi, yi, ti in zip(x_back,y_back,t_back): #loop over vector length
                    if (xi, yi) not in hit_dict_b: #getting the hit times of each fibre
                        hit_dict_b[(xi, yi)] = ti
                    else:
                        hit_dict_b[(xi, yi)] = min(hit_dict_b[(xi, yi)], ti) 

                h_back  = np.zeros(bins)
                for (xpos, ypos), tiso in hit_dict_b.items():
                    xi = np.digitize(xpos, xedges1) - 1
                    yi = np.digitize(ypos, yedges1) - 1
                    h_back[xi, yi] = tiso
                

            elif self.med_time_radio.isChecked():
                print("Median Time selected")

            # Mask the bins with zero entries
            h_back_m = np.ma.masked_where(h_back == 0, h_back)
            h_front_m = np.ma.masked_where(h_front == 0, h_front)

            cmap = plt.cm.viridis
            cmap.set_bad(color='white')  # Set color for masked values

            # Determine the common color scale range
            vmin = min(h_back_m.min(), h_front_m.min())
            vmax = max(h_back_m.max(), h_front_m.max())

            # Plot the 2D histograms
            self.ax1.clear()
            self.ax1.imshow(h_back_m.T, origin='lower', aspect='auto', extent=[xedges1[0], xedges1[-1], yedges1[0], yedges1[-1]], cmap=cmap, norm=LogNorm(vmin=vmin, vmax=vmax))
            self.ax1.set_title('Back panel')

            self.ax2.clear()
            im2 = self.ax2.imshow(h_front_m.T, origin='lower', aspect='auto', extent=[xedges2[0], xedges2[-1], yedges2[0], yedges2[-1]], cmap=cmap, norm=LogNorm(vmin=vmin, vmax=vmax))
            self.ax2.set_title('Front panel')
            if self.colorbar2:
                self.colorbar2.remove()
            self.colorbar2 = self.figure2.colorbar(im2, ax=self.ax2, norm=LogNorm(vmin=vmin, vmax=vmax))

            # Draw dotted circles on the histograms
            circle1 = plt.Circle((0, 0), radius=900, color='red', fill=False, linestyle='dotted')
            circle2 = plt.Circle((0, 0), radius=900, color='red', fill=False, linestyle='dotted')

            self.ax1.add_patch(circle1)
            self.ax2.add_patch(circle2)

            self.figure1.canvas.draw()
            self.figure2.canvas.draw()
        except Exception as e:
            print(f"An error occurred: {e}")

if __name__ == '__main__':
    if len(sys.argv) != 2:
        print("Usage: python script_name.py <input_file>")
        sys.exit(1)
    
    input_file = sys.argv[1]
    app = QApplication([])
    window = EventDisplay(input_file)
    window.show()
    app.exec_()