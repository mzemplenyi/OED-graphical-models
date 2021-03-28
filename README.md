# Description
This repository contains the code to implement the Bayesian optimal experimental design algorithm described in:
Michele Zemplenyi and Jeffrey W. Miller. (2021) "Bayesian optimal experimental design for inferring causal structures graphical models." 

This codebase builds upon the BDAGL (Bayesian Directed Acyclic Graph Learning) Matlab/C/Java package developed by Daniel Eaton and Kevin Murphy, which supports Bayesian inference about DAG structures using dynamic programming and MCMC. Eaton and Murphy distributed their code under the GNU Lesser General Public License. Our code is also distributed under the [GNU Lesser General Public License](https://www.gnu.org/licenses/lgpl-3.0.en.html). See the [BDAGL website](https://www.cs.ubc.ca/~murphyk/Software/BDAGL/) for more information about the BDAGL package. We are grateful to Professor Eaton and Professor Murphy for making their software publicly available.

## Getting started 
All relevant code is included in this repository. The code has only been tested on Matlab R2017a and is not guaranteed to work on other versions. To use this code:
* Clone this repository.
* In Matlab, change the directory to the location of the 'BayesianOED' folder.
* Run the script 'mkPath', which adds all necessary folders to the path and compiles Mex files.
* For simulations:
  * Specify the structure and conditional probability distribution for the DAG that will be used to generate data in the file 'mkBnet'.
  * Enter the desired simulation settings in the 'startSim' file and then run that script.  
