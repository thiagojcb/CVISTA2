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
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar

class EventDisplay(QMainWindow):
    def __init__(self, input_file, evt_no = None):
        super().__init__()
        self.input_file = input_file
        self.evt_no = evt_no
        self.initUI()
        self.colorbar = None
        self.reverse_pmt_dict = {}
        self.pmt_dict = {}
        self.half_length = 0

        # some geometry var. hard coded for now.
        self.scint_radius      = 900
        
        # Plot data if event number is provided
        if self.evt_no is not None:
            self.entry_spinbox.setValue(self.evt_no)
            self.plot_data()

    def initUI(self):
        self.setWindowTitle('Event Display')
        self.setGeometry(100, 100, 1200, 600)  # Adjusted window size

        central_widget = QWidget()
        self.setCentralWidget(central_widget)
        layout = QVBoxLayout(central_widget)

        hBox = QHBoxLayout()

        # MC group box
        group_MC_box = QGroupBox("Truth")
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
        group_Data_box = QGroupBox("Readout")
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
        self.canvas1 = FigureCanvas(self.figure1)
        self.canvas2 = FigureCanvas(self.figure2)
        self.hist_layout.addWidget(self.canvas1)
        self.hist_layout.addWidget(self.canvas2)

        # Add toolbars
        self.toolbar_layout = QHBoxLayout()
        self.toolbar1 = NavigationToolbar(self.canvas1, self)
        self.toolbar2 = NavigationToolbar(self.canvas2, self)
        self.toolbar_layout.addWidget(self.toolbar1)
        self.toolbar_layout.addWidget(self.toolbar2)
        layout.addLayout(self.toolbar_layout)

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

    def on_pick(self, event):
        ind = event.ind[0]
        x = event.artist.get_offsets()[ind, 0]
        y = event.artist.get_offsets()[ind, 1]
        # c = event.artist.get_array()[ind]
        z = float(event.artist.get_gid())
        pmt_id = self.get_pmtID(x,y,z)
        self.ax3.clear()
        times = self.hit_dict.get((x,y,z),None)
        counts, bin_edges, patches = self.ax3.hist(times, range(0,250))
        leg3_text  = f'SiPM#{pmt_id} ({x:.2f},{y:.2f},{z:.2f}) mm\n'
        leg3_text += f'{len(times)} PEs'
        style = dict(size=8, color='gray')
        self.ax3.text(120,1.06*np.max(counts),leg3_text,**style)
        self.ax3.set_xlabel('PE time (ns)')
        self.ax3.set_ylabel('Entries')
        self.figure3.canvas.draw()

    def update_annot(self, ind, scatter, annot):
        pos = scatter.get_offsets()[ind["ind"][0]]
        annot.xy = pos
        text = f"X: {pos[0]:.2f}\nY: {pos[1]:.2f}\nValue: {scatter.get_array()[ind['ind'][0]]:.2f}"
        annot.set_text(text)
        annot.get_bbox_patch().set_alpha(0.4)

    def hover(self, event):
        vis1 = self.annot1.get_visible()
        vis2 = self.annot2.get_visible()
        if event.inaxes == self.ax1:
            scatter = self.ax1.collections[0]
            cont, ind = scatter.contains(event)
            if cont:
                self.update_annot(ind, scatter, self.annot1)
                self.annot1.set_visible(True)
                self.figure1.canvas.draw_idle()
                return
        elif event.inaxes == self.ax2:
            scatter = self.ax2.collections[0]
            cont, ind = scatter.contains(event)
            if cont:
                self.update_annot(ind, scatter, self.annot2)
                self.annot2.set_visible(True)
                self.figure2.canvas.draw_idle()
                return
        if vis1:
            self.annot1.set_visible(False)
            self.figure1.canvas.draw_idle()
        if vis2:
            self.annot2.set_visible(False)
            self.figure2.canvas.draw_idle()

    def get_pmtID(self,x,y,z):
        coords = (x, y, z)
        return self.reverse_pmt_dict.get(coords, None)

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

            if z[0] < 0:
                self.half_length = -1*z[0]
            else:
                self.half_length = z[0]

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
            self.hit_dict = {}
            for xi, yi, zi, ti in zip(x,y,z,t): #loop over vector length
                if (xi, yi, zi) not in self.hit_dict: #getting the hit times of each fibre
                    self.hit_dict[(xi, yi, zi)] = [ti]
                else:
                    self.hit_dict[(xi, yi, zi)].append(ti)

            # Create a dictionary to map pmtid to their corresponding values
            self.pmt_dict = {id: (x, y, z) for id, x, y, z in zip(pmtID, pmtX, pmtY, pmtZ)}
            self.reverse_pmt_dict = {coords: id for id, coords in self.pmt_dict.items()}

            if self.signal_radio.isChecked():

                # Extract the corresponding pmtx, pmty, and pmtz values for the flagged pmtids
                flagged_x = np.array([self.pmt_dict[id][0] for id in mcPMTID])
                flagged_y = np.array([self.pmt_dict[id][1] for id in mcPMTID])
                flagged_z = np.array([self.pmt_dict[id][2] for id in mcPMTID])

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

                # Get times with values from the dictionary
                front_npe_x = np.array([])
                front_npe_y = np.array([])
                front_npe   = np.array([])
                back_npe_x = np.array([])
                back_npe_y = np.array([])
                back_npe   = np.array([])
                for (xi, yi, zi), tiso in self.hit_dict.items():
                    if zi > 0:
                        front_npe_x = np.append(front_npe_x,xi)
                        front_npe_y = np.append(front_npe_y,yi)
                        front_npe   = np.append(front_npe,min(tiso))
                    else:
                        back_npe_x = np.append(back_npe_x,xi)
                        back_npe_y = np.append(back_npe_y,yi)
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
            sizes = 5
            if self.time_radio.isChecked():
                scatter1 = self.ax1.scatter(back_npe_x, back_npe_y, c=back_npe, s=sizes, alpha=0.5, norm=norm, picker=True)
            else:
                scatter1 = self.ax1.scatter(back_npe_x, back_npe_y, c=back_npe, s=sizes, alpha=0.5, vmin=vmin, vmax=vmax, picker=True)
            scatter1.set_gid(str(-1*self.half_length)) # HARD CODE Z!!!
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
            # sizes = 5
            if self.time_radio.isChecked():
                scatter2 = self.ax2.scatter(front_npe_x, front_npe_y, c=front_npe, s=sizes, alpha=0.5, norm=norm, picker=True)
            else:
                scatter2 = self.ax2.scatter(front_npe_x, front_npe_y, c=front_npe, s=sizes, alpha=0.5, vmin=vmin, vmax=vmax, picker=True)
            scatter2.set_gid(str(self.half_length)) # HARD CODE Z!!!
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
            circle1 = plt.Circle((0, 0), radius=self.scint_radius, color='red', fill=False, linestyle='dotted')
            circle2 = plt.Circle((0, 0), radius=self.scint_radius, color='red', fill=False, linestyle='dotted')

            self.ax1.add_patch(circle1)
            self.ax2.add_patch(circle2)

            self.figure1.canvas.draw()
            self.figure2.canvas.draw()

            # Create an annotation object
            self.annot1 = self.ax1.annotate("", xy=(0,0), xytext=(20,20),
                                           textcoords="offset points",
                                           bbox=dict(boxstyle="round", fc="w"),
                                           arrowprops=dict(arrowstyle="->"))
            self.annot1.set_visible(False)

            self.annot2 = self.ax2.annotate("", xy=(0,0), xytext=(20,20),
                                           textcoords="offset points",
                                           bbox=dict(boxstyle="round", fc="w"),
                                           arrowprops=dict(arrowstyle="->"))
            self.annot2.set_visible(False)

            # Connect the motion_notify_event to update the annotation
            self.figure1.canvas.mpl_connect("motion_notify_event", self.hover)
            self.figure2.canvas.mpl_connect("motion_notify_event", self.hover)

            # Connect the pick event to a handler function
            self.figure1.canvas.mpl_connect('pick_event', self.on_pick)
            self.figure2.canvas.mpl_connect('pick_event', self.on_pick)

        except Exception as e:
            print(f"An error occurred: {e}")

if __name__ == '__main__':
    if len(sys.argv) < 2 or len(sys.argv) > 3:
        print("Usage: python script_name.py <input_file> [evt_no]")
        sys.exit(1)
    
    input_file = sys.argv[1]
    evt_no = int(sys.argv[2]) if len(sys.argv) == 3 else None
    app = QApplication([])
    window = EventDisplay(input_file, evt_no)
    window.show()
    app.exec_()