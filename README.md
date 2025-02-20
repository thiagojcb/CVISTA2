# CVISTA2
A data visualisation GUI for **ntuple** files of opaque scintillator detectors simulation based on [Ratpac-two](https://github.com/rat-pac/ratpac-two) simulation and analysis package.

## Dependencies
- uproot
- numpy
- awkward
- matplotlib
- particle
- PyQt5
  
Use `pip install <module>` if needed.

## Usage
`python CVISTA2.py <input_file> [event/entry number]`
- `input_file` is a ntuple output file of Ratpac.
- `event number` is an optional argument. If given, the GUI will display this event at start up.
