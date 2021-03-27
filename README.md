# OED-graphical-models
This repository contains the code to implement the Bayesian optimal experimental design algorithm described in:
Michele Zemplenyi and Jeffrey W. Miller. (2021) "Bayesian optimal experimental design for inferring causal structures graphical models." 

This codebase builds upon the BDAGL (Bayesian Directed Acyclic Graph Learning) Matlab/C/Java package developed by Daniel Eaton and Kevin Murphy, which supports Bayesian inference about DAG structures using dynamic programming and MCMC. Eaton and Murphy distributed their code under the Lesser (formerly Library) GNU Public License. See here for more information about the BDAGL package: https://www.cs.ubc.ca/~murphyk/Software/BDAGL/. We are grateful to Professors Eaton and Murphy for making their software publicly available.

All relevant code is included in this repository. The code has only been tested on Matlab R2017a and is not guaranteed to work on other versions. To use this code:
-Clone this repository or download all the code.
-In Matlab, change the directory to the location of the 'Bayesian_OED' folder.
-Run the script 'mkPath', which adds all necessary folders to the path and compiles Mex files. 



