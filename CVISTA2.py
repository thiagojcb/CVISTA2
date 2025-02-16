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
from PyQt5.QtWidgets import QHBoxLayout, QSpinBox, QLabel, QRadioButton, QGroupBox, QTextEdit
from matplotlib.gridspec import GridSpec
from mpl_toolkits.axes_grid1 import make_axes_locatable


class EventDisplay(QMainWindow):
    def __init__(self, input_file):
        super().__init__()
        self.input_file = input_file
        self.initUI()
        self.colorbar = None

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

        # Create a layout for the SiPMs
        self.hist_layout = QHBoxLayout()
        layout.addLayout(self.hist_layout)

        # Create a figure with shared y-axis
        self.figure1, self.ax1 = plt.subplots()
        self.figure2, self.ax2 = plt.subplots()

        # Add figures to the layout
        self.hist_layout.addWidget(self.figure1.canvas)
        self.hist_layout.addWidget(self.figure2.canvas)

        control_layout = QHBoxLayout()
        button_layout  = QVBoxLayout()
        display_layout = QHBoxLayout()

        # Event selector
        self.entry_label = QLabel('Select Entry:')
        display_layout.addWidget(self.entry_label)

        self.entry_spinbox = QSpinBox()
        self.entry_spinbox.setMinimum(0)
        display_layout.addWidget(self.entry_spinbox)

        # Display button
        self.plot_button = QPushButton('Display event')
        self.plot_button.clicked.connect(self.plot_data)
        display_layout.addWidget(self.plot_button)

        button_layout.addLayout(display_layout)

        # Next button
        self.next_button = QPushButton('Display next')
        self.next_button.clicked.connect(self.display_next_event)
        display_layout.addWidget(self.next_button)
        
        # Previous button
        self.prev_button = QPushButton('Display previous')
        self.prev_button.clicked.connect(self.display_prev_event)
        display_layout.addWidget(self.prev_button)

        # Text with basic event info
        self.text_edit = QTextEdit()
        button_layout.addWidget(self.text_edit)

        control_layout.addLayout(button_layout)
        
        self.figure3, self.ax3 = plt.subplots()
        control_layout.addWidget(self.figure3.canvas)
        
        layout.addLayout(control_layout)

    def display_next_event(self):
        current_index = self.entry_spinbox.value()
        if current_index < self.entry_spinbox.maximum():
            self.entry_spinbox.setValue(current_index + 1)
            self.plot_data()
            
    def display_prev_event(self):
        current_index = self.entry_spinbox.value()
        if current_index > self.entry_spinbox.minimum():
            self.entry_spinbox.setValue(current_index - 1)
            self.plot_data()

    def plot_data(self):
        try:
            # Read data from ROOT file
            file = uproot.open(self.input_file)
            tree = file['output']
            meta = file['meta']
            entry_index = self.entry_spinbox.value()
            nentries = tree.num_entries
            self.entry_spinbox.setMaximum(nentries-1)

            # check if branches exist
            branches = tree.keys()
            reco_exists = 'x_FitCentroid' in branches

            pmtID = meta["pmtId"].array()[0]
            pmtID = ak.to_numpy(pmtID)
            pmtX  = meta["pmtX"].array()[0]
            pmtX = ak.to_numpy(pmtX)
            pmtY  = meta["pmtY"].array()[0]
            pmtY = ak.to_numpy(pmtY)
            pmtZ  = meta["pmtZ"].array()[0]
            pmtZ = ak.to_numpy(pmtZ)
     
            # Initialize text
            summary_text = 'File: '+self.input_file+f' ({nentries} entries)\n'
            summary_text += f'Entry: {entry_index}\n'

            mcpdg    = tree['mcpdg'].array(entry_start=entry_index, entry_stop=entry_index+1)[0]
            mcpdg    = ak.to_numpy(mcpdg)
            summary_text += f"MC particle: {mcpdg[0]}\n"

            mcke     = tree['mcke'].array(entry_start=entry_index, entry_stop=entry_index+1)[0]
            mcke     = ak.to_numpy(mcke)
            summary_text += f"KE: {mcke[0]:.1f} MeV\n"

            x_th = tree['mcx'].array(entry_start=entry_index, entry_stop=entry_index+1)[0]
            x_th = ak.to_numpy(x_th)
            summary_text +=  f"X: {x_th[0]:.1f} mm, "
            y_th = tree['mcy'].array(entry_start=entry_index, entry_stop=entry_index+1)[0]
            y_th = ak.to_numpy(y_th)
            summary_text += f"Y: {y_th[0]:.1f} mm, "
            z_th = tree['mcz'].array(entry_start=entry_index, entry_stop=entry_index+1)[0]
            z_th = ak.to_numpy(z_th)
            summary_text += f"Z: {z_th[0]:.1f} mm\n"

            x = tree['mcPEx'].array(entry_start=entry_index, entry_stop=entry_index+1)[0]
            x = ak.to_numpy(x)

            y = tree['mcPEy'].array(entry_start=entry_index, entry_stop=entry_index+1)[0]
            y = ak.to_numpy(y)

            z = tree['mcPEz'].array(entry_start=entry_index, entry_stop=entry_index+1)[0]
            z = ak.to_numpy(z)

            t = tree['mcPEHitTime'].array(entry_start=entry_index, entry_stop=entry_index+1)[0]
            t = ak.to_numpy(t)

            # mcpecount= tree['mcpecount'].array(entry_start=entry_index, entry_stop=entry_index+1)[0]
            # mcpecount= ak.to_numpy(mcpecount)

            # mcnhits  = tree['mcnhits'].array(entry_start=entry_index, entry_stop=entry_index+1)[0]
            # mcnhits  = ak.to_numpy(mcnhits)

            mcPMTID = tree['mcPMTID'].array(entry_start=entry_index, entry_stop=entry_index+1)[0]
            mcPMTID = ak.to_numpy(mcPMTID)

            mcPMTNPE = tree['mcPMTNPE'].array(entry_start=entry_index, entry_stop=entry_index+1)[0]
            mcPMTNPE = ak.to_numpy(mcPMTNPE)

            mcpecount = sum(mcPMTNPE)
            mcnhits   = len(mcPMTNPE>0)

            summary_text += f"SiPMs w/ PEs: {mcnhits},"
            summary_text += f"Total PEs: {mcpecount}\n"

            if reco_exists:
                x_reco = tree['x_FitCentroid'].array(entry_start=entry_index, entry_stop=entry_index+1)[0]
                x_reco = ak.to_numpy(x_reco)
                summary_text += f"X reco: {x_reco[0]:.1f} mm, "
                y_reco = tree['y_FitCentroid'].array(entry_start=entry_index, entry_stop=entry_index+1)[0]
                y_reco = ak.to_numpy(y_reco)
                summary_text += f"Y reco: {y_reco[0]:.1f} mm, "
                z_reco = tree['z_FitCentroid'].array(entry_start=entry_index, entry_stop=entry_index+1)[0]
                z_reco = ak.to_numpy(z_reco)
                summary_text += f"Z reco: {z_reco[0]:.1f} mm\n"
            else:
                summary_text += "X reco: NaN, "
                summary_text += "Y reco: NaN, "
                summary_text += "Z reco: NaN\n"

            self.text_edit.setText(summary_text)

            # Define fixed ranges for the histograms
            x_range = (-1250, 1250)
            y_range = (-1250, 1250)

            if self.signal_radio.isChecked():
                # Create a dictionary to map pmtid to their corresponding values
                pmt_dict = {id: (x, y, z) for id, x, y, z in zip(pmtID, pmtX, pmtY, pmtZ)}

                # Extract the corresponding pmtx, pmty, and pmtz values for the flagged pmtids
                flagged_x = np.array([pmt_dict[id][0] for id in mcPMTID])
                flagged_y = np.array([pmt_dict[id][1] for id in mcPMTID])
                flagged_z = np.array([pmt_dict[id][2] for id in mcPMTID])

                back_npe   = mcPMTNPE[flagged_z < 0]
                back_npe_x = flagged_x[flagged_z < 0]
                back_npe_y = flagged_y[flagged_z < 0]
                withPE    = back_npe>0

                back_npe   = back_npe[withPE]
                back_npe_x = back_npe_x[withPE]
                back_npe_y = back_npe_y[withPE]

                front_npe   = mcPMTNPE[flagged_z > 0]
                front_npe_x = flagged_x[flagged_z > 0]
                front_npe_y = flagged_y[flagged_z > 0]
                withPE    = front_npe>0

                front_npe   = front_npe[withPE]
                front_npe_x = front_npe_x[withPE]
                front_npe_y = front_npe_y[withPE]

            elif self.time_radio.isChecked():
                # Time MC data plotting
                # get back and front channels
                negative_pez_indices = z < 0
                positive_pez_indices = z >= 0

                x_back = x[negative_pez_indices]
                y_back = y[negative_pez_indices]
                t_back = t[negative_pez_indices]

                x_front = x[positive_pez_indices]
                y_front = y[positive_pez_indices]
                t_front = t[positive_pez_indices]

                hit_dict_f = {}
                for xi, yi, ti in zip(x_front,y_front,t_front): #loop over vector length
                    if (xi, yi) not in hit_dict_f: #getting the hit times of each fibre
                        hit_dict_f[(xi, yi)] = [ti]
                    else:
                        hit_dict_f[(xi, yi)].append(ti)

                # Get times with values from the dictionary
                front_npe_x = np.array([])
                front_npe_y = np.array([])
                front_npe   = np.array([])
                for (xpos, ypos), tiso in hit_dict_f.items():
                    front_npe_x = np.append(front_npe_x,xpos)
                    front_npe_y = np.append(front_npe_y,ypos)
                    front_npe   = np.append(front_npe,min(tiso))
                
                hit_dict_b = {}
                for xi, yi, ti in zip(x_back,y_back,t_back): #loop over vector length
                    if (xi, yi) not in hit_dict_b: #getting the hit times of each fibre
                        hit_dict_b[(xi, yi)] = [ti]
                    else:
                        hit_dict_b[(xi, yi)].append(ti) 

                back_npe_x = np.array([])
                back_npe_y = np.array([])
                back_npe   = np.array([])
                for (xpos, ypos), tiso in hit_dict_b.items():
                    back_npe_x = np.append(back_npe_x,xpos)
                    back_npe_y = np.append(back_npe_y,ypos)
                    back_npe   = np.append(back_npe,min(tiso))

            elif self.med_time_radio.isChecked():
                print("Median Time selected")

            cmap = plt.cm.viridis
            cmap.set_bad(color='white')  # Set color for masked values

            # Determine the common color scale range
            vmin = min(back_npe.min(), front_npe.min())
            vmax = max(back_npe.max(), front_npe.max())

            norm = LogNorm(vmin=vmin, vmax=vmax)

            # Plot the scatter 
            self.ax1.clear()
            # sizes = back_npe
            sizes = 1
            if self.time_radio.isChecked():
                scatter1 = self.ax1.scatter(back_npe_x, back_npe_y, c=back_npe, s=sizes, alpha=0.5, norm=norm)
            else:
                scatter1 = self.ax1.scatter(back_npe_x, back_npe_y, c=back_npe, s=sizes, alpha=0.5, vmin=vmin, vmax=vmax)
            self.ax1.set_xlim(x_range)
            self.ax1.set_ylim(y_range)

            self.ax1.set_xlabel('X (mm)')
            self.ax1.set_ylabel('Y (mm)')
            leg1_text  = '-Z Channels\n'
            leg1_text += f'{len(back_npe_x)} SiPMs\n'
            if self.time_radio.isChecked():
                leg1_text += f'{min(back_npe):.1f} ns (1st PE)'
            else:
                leg1_text += f'{sum(back_npe)} PEs'
            style = dict(size=8, color='gray')
            self.ax1.text(800,600,leg1_text,**style)

            self.ax2.clear()
            # sizes = front_npe
            sizes = 1
            if self.time_radio.isChecked():
                scatter2 = self.ax2.scatter(front_npe_x, front_npe_y, c=front_npe, s=sizes, alpha=0.5, norm=norm)
            else:
                scatter2 = self.ax2.scatter(front_npe_x, front_npe_y, c=front_npe, s=sizes, alpha=0.5, vmin=vmin, vmax=vmax)
            self.ax2.set_xlim(x_range)
            self.ax2.set_ylim(y_range)
            self.ax2.set_xlabel('X (mm)')
            leg2_text  = '+Z Channels\n'
            leg2_text += f'{len(front_npe_x)} SiPMs\n'
            if self.time_radio.isChecked():
                leg2_text += f'{min(front_npe):.1f} ns (1st PE)'
            else:
                leg2_text += f'{sum(front_npe)} PEs'
            self.ax2.text(800,600,leg2_text,**style)

            if self.colorbar:
                self.colorbar.remove()
            divider = make_axes_locatable(self.ax2)
            cax = divider.append_axes("right", size="5%", pad=0.05)
            if self.time_radio.isChecked():
                self.colorbar = self.figure2.colorbar(scatter2, cax=cax, norm=norm)
            else:
                # self.colorbar = self.figure.colorbar(im2, cax=cax)
                self.colorbar = self.figure2.colorbar(scatter2, cax=cax)

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