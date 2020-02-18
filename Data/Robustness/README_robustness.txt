Network_<gamma>.mat (2/17/19): contains GLASSO networks for the robustness subset (14 patients both ictal and preictal). <gamma> is the value of the EBIC parameter
    - wSim_<beta>: similarity newtwork with distance weighting parameter beta


Partitions_<gamma>.mat(2/17/19): contains partition info for the robustness subset (12 patients both ictal and preictal) at multiple beta values.
 <gamma> is the value of the EBIC parameter gamma used in generating the GLASSO networks. 
    - beta: distance weighting parameter beta applied before modularity maximization
    NOTE: Partitions have only been calculated at beta=0.01 for patients Study023 and HUP106 (ictal and preictal)
