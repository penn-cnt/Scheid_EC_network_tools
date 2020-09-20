# EC_network_tools
Tools for effective connectivity network generation and analysis

## Requirements
- Matlab >= 2018a
- R >= 3.0.0
- R libraries: bootnet, R.matlab
- The commandline tool "Rscript" must be installed and should be on the system PATH
- GenLouvain Version 2.2: https://github.com/GenLouvain/GenLouvain

Note, code in this repository has only been tested on linux and mac systems and may not be compatible with other environments. 

## Example

```matlab 
addpath(helper_functions)

% Get regularized partial correlation matrices across windowed data
timeseries_data = [48 channels x 5000 samples]
windowSize = 500;
Networks = getNets('./', timeseries_data, windowSize, 0.5, [0, 0.1, .5])

% Partition Network similarity matrix into desired number of communities
init_gamma = (0.8:0.05:1); % Set of initial spatial resolution parameters 
nTarget = 3;               % Target number of communities 
Qiter = 100;               % iterations of Louvain algorithm
Partitions = findSimComms(Networks.wSim_0_5, init_gamma, Qiter, nTarget)

```
