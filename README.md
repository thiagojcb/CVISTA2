# CVISTA2
A data visualisation GUI for **ntuple** files of opaque scintillator detectors simulation based on [Ratpac-two](https://github.com/rat-pac/ratpac-two) simulation and analysis package.
![image](https://github.com/user-attachments/assets/948e44fa-da53-40c9-a0f6-c82fb5ffff1d)

## Dependencies
- uproot
- numpy
- awkward
- matplotlib
- particle
- PyQt5
  
Use `pip install <module>` if needed.

Developed and tested with Python 3.11.6.

## Usage
`python CVISTA2.py <input_file> [event/entry number]`
- `input_file` is a ntuple output file of Ratpac.
- `event number` is an optional argument. If given, the GUI will display this event at start up.

### Radio buttons
The radio buttons on the top of the GUI indicates what kind of information you want to be plotted on the histograms.
- nPE: the number of simulation PEs collected by each SiPM.
- First Hit Time: the simulated time of the first PE to reach a SiPM.
- SiPM Charge: the integrated charge of each SiPM
- SiPM Time: the pulse start time of each SiPM
