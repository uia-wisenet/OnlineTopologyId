# OnlineTopologyId
Codes for replicating results in paper "Online Topology Identification from Vector Autoregressive Time Series" by Zaman, Lopez-Ramos, Romero and Beferull-Lozano

# Introduction

The folder gsim_onlineTopID contains the MATLAB code for generating the figures of the simulation results in the paper

B. Zaman, L. M. Lopez-Ramos, D. Romero, B. Beferull-Lozano, “Online Topology Identification from Vector Autoregressive Time-series”, To appear in IEEE Transactions on Signal Processing. 


## Guidlines
To run the code the experiments, please follow the following steps after downloading the code:
* enter the folder **gsim_onlineTopID** and run ```gsimStartup```. This will generate ```gsim.m```in **gsim_onlineTopID**.
* The file ```Experiments\OnlineTopIdTSPExperiments.m``` contains all the experiments in the paper.
* The first digit of the epxeriment number corresponds to the Figure number in the paper.
* One can see the saved experiment results in **Experiments\OnlineTopIdTSPExperiments** or can execute the experiment. For instance, if one wants to see the already saved Fig. 2 only, run ``` gsim(1,20)```. To run the experiment, execute ```gsim(0,20,[])```.

## Citation
```
@article{zaman2020online,
  title={Online Topology Identification from Vector Autoregressive Time-series},
  author={B. Zaman and L. M. Lopez-Ramos and D. Romero and B. Beferull-Lozano},
  journal={to appear in IEEE Transactions on Signal Processing},
  }
