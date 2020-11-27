# OnlineTopologyId
Codes for replicating results in paper:

B. Zaman, L. M. Lopez-Ramos, D. Romero, B. Beferull-Lozano, “Online Topology Identification from Vector Autoregressive Time-series”, To appear in IEEE Transactions on Signal Processing. 

## Introduction
The folder gsim_onlineTopID contains the MATLAB code for generating the figures of the simulation results in the paper.


## Guidelines
To run the code to replicate the experiments, please follow these steps after downloading the code:
* enter the folder **gsim_onlineTopID** and run ```gsimStartup```. This will generate ```gsim.m```in **gsim_onlineTopID**.
* The file ```Experiments\OnlineTopIdTSPExperiments.m``` contains all the experiments in the paper. Each experiment is a method identified with a number corresponding to the figure it generates.

* One can execute the experiment, or see the saved experiment results (without actually re-running the experiment) in **Experiments\OnlineTopIdTSPExperiments**. The first argument passed to gsim is an "onlyplot" flag:
  * To see the already saved Fig. 2 only, run ``` gsim(1, 2)```. 
  * To run the experiment, execute ```gsim(0, 2)```.

## Citation
```
@article{zaman2020online,
  title={Online Topology Identification from Vector Autoregressive Time-series},
  author={B. Zaman and L. M. Lopez-Ramos and D. Romero and B. Beferull-Lozano},
  journal={to appear in IEEE Transactions on Signal Processing},
  }
