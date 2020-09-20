# EC_network_tools
Tools for effective connectivity network generation and analysis

## Requirements
- Matlab >= 2018a
- R >= 3.0.0
- R libraries: bootnet, R.matlab
- The commandline tool "Rscript" must be installed and should be on the system PATH

Note, code in this repository has only been tested on linux and mac systems and may not be compatible with other environments. 

## Example

```matlab 
addpath(helper_functions)

timeseries_data = [48 channels x 5000 samples]
windowSize = 500;

Networks = getNets('./', timeseries_data, windowSize, 0.5, [0, 0.1, .5])
```
