# 4DNucleomeEquation

4DNucleomeEquation simulation analysis and statistical inference of gene expression regulated by enhancer-promoter interaction across space and time described in the paper "4D nucleome equation predicts gene expression controlled by long-range enhancer-promoter interaction".

## **Environment setup**

First let’s set up the environment to run the simulations. The simulated system is complex and requires parallel computing. We can put the simulated code on the server for computing, or test a small demo on the local machine. First of all, download the repository as a zip file clicking the green button in the top right of this page to your local machine.

### **Download the scripts**

First of all, download the repository as a zip file clicking the green button in the top right of this page to your local machine.

### **Requirements**

The software prerequisites for the package to work are Matlab R2018a or later version of Matlab. You can start parallel pool (parpool) using the command `parpool`. Then, you can obtain the number of workers.

```
%% current parallel pool
if isempty(gcp('nocreate'))
    parpool(4); % 4 is the number of workers.It can be adjusted according to different local devices and servers
end
```

### **Loading the scripts**

The next step is to manually add the scripts in the same folder.

## Directories

There are two folders, and each folder contains several .m files. Next, we give detailed explanation each script:

### SimulationAnalysis

####  testdG.m, testKEP.m

The script can simulate and analyze different values of enhancer-promoter (E-P) genome distance (or E-P interaction strength). We simulate multiple times in parallel and store the simulated data in the folder generated according to the parameters.

#### ParametersBurst.m

This function generates a parameter structure to be passed to the simulation framework. The parameters are saved fields of a structure called `params`. We need to pass in the user defined parameters and set some default parameter values.

#### SimulateBurst.m

This .m files is core simulation framework. This script includes pre-allocating variable, initializing model, updating model and saving results.

#### InitializeConnectivityMatrix.m

This function is used to initialize generalized Rouse model, and was designed for initialized params structure input format. We supply this file to generate connectivity matrix to represent the connection of linear monomers. This file should not be used for other purposes.

#### AnalyseBurst.m

This .m file is the core analysis process. It contains two parts: Statistical analysis of simulation data and theoretical calculation of statistical indicators. We can obtain the results numerically and theoretically containing: mean mRNA level, CV, etc.

#### Poissbeta.m

This .m file can theoretically calculate the probability density function of ON-OFF model.

#### drawMeanVardG.m, drawMeanVarKEP.m

This script can plot the mean and CV changes of E-P genome distance (or E-P interaction strength) based on simulation results. 

#### drawBC.m

This script theoretically computes the bimodal coefficient and determine the peak number boundaries.

#### drawDGpeakinformation.m, drawKEPpeakinformation.m

This script theoretically calculate the peak information for E-P genome distance (E-P interaction strength) under a fixed E-P interaction strength (E-P genome distance). We show the evolutionary process of the peak numbers and peak probabilities.

### StatisticalInference

#### dataToFit.mat, dGvsEncounterProb.mat

`dataToFit.mat`  provides the experimental data we are going to fit. The data includes the E-P genome distance, mRNA distribution, E-P encounter prob and mRNA mean level for different cell lines.  `dGvsEncounterProb.mat`  provides the encounter prob and corresponding E-P genomic distance. 

#### fitmRNADistribution.m

This script  the master program for estimating the  E-P interaction dynamics and gene expression dynamics parameters from smFISH data and E-P genomic distance data.

#### fitEPEncounterProb.m

This script compares theoretical results based on statistically inferred parameters with experimental data. We show the relationship between E-P genome distance and contact probability. Mean eGFP mRNA level plotted against contact probability between the ectopic Sox2 promoter and SCR insertions. CV of eGFP level against contact probability.

#### CEInferenceBurst.m

This script is the core code of statical inference. We use the minimum cross entropy method to estimate the parameters.

#### CalculateTheoryProb.m

This .m file theoretical calculate of mRNA distribution.
