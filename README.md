### README
This repository contains sample MATLAB code for carrying out global sensitivity analysis to identify key engineering interventions for designing co-culture dynamics, as detailed in our paper
**A molecular toolkit of cross-feeding strains for engineering synthetic yeast communities** *to appear in Nature Microbiology*
by Huadong Peng, Alexander P. S. Darlington, Eric J. South, Hao-Hong Chen, Wei Jiang and Rodrigo Ledesma Amaro.
Correspondence should be addressed to Rodrigo Ledesma-Amaro (r.ledesma-amaro (AT) imperial.ac.uk). 
The repository was written by Alexander P. S. Darlington (a.darlington.1 (AT) warwick.ac.uk). 

### FILE LIST
The repository contains the following files:

##### Sample results
**FIGURE_ode1_and_ode2_strains_design.m** simulates the dynamics of the one and two strain models for the nominal parameter set.
**EFAST_ode2strain** --- folder containing exemplar results of GSA for the two strain coculture model.
**EFAST_ode2MEv1strain** --- folder containing exemplar results of GSA for the two strain coculture augmented with the metabolic pathway model used to inform the DoL experiments.
**EFAST_ode2MEv1strain_vary_ME** --- folder containing exemplar results of GSA for the two strain coculture augmented with the metabolic pathway model used to inform the DoL experiments with the metabolic pathway parameters also varied.
**EFAST_ode2TOXstrain** --- folder containing exemplar results of GSA for the two strain coculture model with metabolite toxicity.
**EFAST_ode2REUPstrain** --- folder containing exemplar results of GSA for the two strain coculture model with metabolite re-uptake.
**EFAST_ode3v1astrain** --- folder containing exemplar results of the GSA for the three strain coculture model with uni-directional metabolite exchange between the strains. 
**EFAST_ode3v1bstrain** --- folder containing exemplar results of the GSA for the three strain coculture model with additional communication between the strains.
**EFAST_ode3v2strain** --- folder containing exemplar results of the GSA for the three strain coculture model with bi-directional metabolite exchange between the strains. 
**EFAST_ode3v2REUPstrain** --- folder containing exemplar results of the GSA for the three strain coculture model with bi-directional metabolite exchange between the strains and metabolite reuptake.         

Sample result figures as SVG files are also provided.

##### Script list
**RUN_efast_#####.m** --- These scripts set up the eFAST sampling problem (including parameter ranges, ODE model choice etc) and carry out the parameter sampling and computationally intensive simulations. The results as saved as *efast_results.mat* in *fname*. These results can then be analysed using **RUN_ANALYSIS.m**.

**RUN_ANALYSIS.m** --- This script loads the sampled parameters and performances from *efast_results.mat* from the folder specificed by *fname* and calcultes the sensitivities for each parameter. It also carries out the T-tests to establish parameters whose sensitivity is more significant than *delta* ("dummy parameter"). The results are saved in the folder *fname* as *efast_stats.mat*. If the user would like to change the statistical significance thresholds this can be achieved by modifying *alpha*. We recommend using Bonferroni correction. Then number of parameters sampled is *length(knames)*.

**FIGURE_#####.m** --- These scripts loads the results from *efast_stats.mat* in the folder specified by *fname* and plot the eFAST sensitivity. Those beginning with "FIGURE_SM_" make the figures available in the supplementary material.

##### Function list
**ode#####strains.m** contains the ODEs for the specified model.
**efast#####strainV2b.m** calls the ODE function, simulates the model for a given parameter set and returns the performance metrics.
**efastRunAnalysisWithConstraintsV2b.m** carries out the parameter sampling and simulations of these parameters needed for 
eFAST analysis.
**efastParameterDistWithConstraints.m** creates the sampled parameters for the model simulations.
**efastSD.m** calculates the sensitivities S_i and S_Ti of each performance metric for each paraemter.
**efastTTest.m** calculates the T-statistics for the S_i and S_Ti values.
**efastSetFreq.m** generates the selection of frequencies needed for the parameter sampling.
**backupcode.m** copies the specified script and needed files to the target analysis folder *fname* as *.txt* files. 

##### SOURCE CODE FOR EFAST ALGORIMTH
For the global sensitivity analysis, we employed the eFAST method developed by Simeone Marino, Ian B Hogue, Christian J Ray and Denise E Kirschner.

The method is fully detailed in their paper **Marino et al. "A methodology for performing global uncertainty and sensitivity analysis in systems biology", Journal of Theoretical Biology, vol. 254, issue. 1., pages 178-196. 2008. DOI [10.1016/j.jtbi.2008.04.011](https://doi.org/10.1016/j.jtbi.2008.04.011).** Their source code is available online at [http://malthus.micro.med.umich.edu/lab/usanalysis.html](http://malthus.micro.med.umich.edu/lab/usanalysis.html).

The functions **efastTTest.m**, **efastSD.m**, **efastSetFreq.m** are renamed but otherwise unchanged and were originally written by Marino et al and retrieved as **efast_ttest.m**, **efast_sd.m**, **SETFREQ.m** available from [their repository](http://malthus.micro.med.umich.edu/lab/usadata/)

The function **efastParameterDistWithConstraints.m** was modified from **parameterdist.m** which was retrieved from [their repository](http://malthus.micro.med.umich.edu/lab/usadata/). The function is modified by Alexander Darlington so that parameters in *X* with indexes *constrainedkidx* sum to 1. This modification enabled us to vary the initial composition of the cocultures without also increasing the initial population size.

The function **efastRunAnalysisWithConstraintsV2b.m** carried out the eFAST analysis for a given model and parameter ranges. This was written by Alexander Darlington using the script **Model_efast.m** originally written by Marino et al. and available from [their repository](http://malthus.micro.med.umich.edu/lab/usadata/).
