#### tax-evasion
# Code reposity for "Tax Evasion and Capital Taxation" by Shahar Rotberg and Joseph Steinberg #

# Data #
Most of our assigned parameter values and calibration targets are taken directly from the literature. There are a few, however, that we constructed ourselves from several data sources.  
<br/>
**Labor productivity process.** We model our idiosyncratic labor productivity process after GKKOC. Their process has two parts: a fixed effect that is constant over an individual's life and transmitted according to an AR(1) process in logs across generations; and a persistent shock that evolves according to an AR(1) process in logs over an individual's life but is not transmitted across generations. Our process is simpler: we have one component that follows one AR(1) process over an individual's life and follows a different AR(1) process across generations. To set the parameters of our process, we simulate data from GKKOC's process and estimate the parameters of our process on that data. Specifically, we first simulate the two GKKOC processes. We then estimate our own process on this simulated data.

**Entrepreneurial opportunity process.** We use the 2016 wave of the Survey of Consumer Finances (SCF) to calculate (i) the fraction of all working-age households with strictly positive business or farm income, and (ii) the fraction of households at age 25 with strictly positive business or farm income. We set the probability of having an entrepreneurial opportunity at birth to (i). Given our assignment of the probability of losing an entrepreneurial opportunity, we then set the probability of regaining an entrepreneurial opportunity so that the model is consistent with (ii).

**Labor income tax rates.** We use data from Table B-1 in https://www.cbo.gov/publication/57404 to compute the average labor income tax rates for different segments of the distribution, and then map these rates to our discretized labor productivity process, which has five states. The bottom 20\% of households, which corresponds to the first two labor productivity states in our model, pay no taxes. The middle 60%, which corresponds to the second state in our model, have a 12.5% average tax rate. The top 20\, which corresponds to the fourth state in our model, pay 24.4% in taxes on average. To compute the tax rate for the fifth state in our model, which contains the top 2.5% of labor income earners, we take the weighted average of the tax rates for the top 1% (30.2%) and the 96th-99th percentiles (24.2%).

**Corporate capital share.** We calculate the corporate income share by dividing profits of nonfinancial domestic corporations (line 13 in NIPA Table 6.16D) by nominal gross domestic product (line 1 of NIPA Table 1.1.5). We use the average corporate income share for 2010--2017 in our calibration to maintain consistency with the timing of the BLS labor share series.

**Detection rate.** To obtain an estimate for the rate at which offshore evasion is detected, we multiply the audit rate for individuals in the top 0.1\% of the income distribution by the probability that an auditor correctly detects an unfulfilled requirement to report offshore wealth. Based on figure A7 in Guyton et al. (2021), we calculate that the average audit rate for individuals in the top 0.1% from 2009-2019 was 8.1%. Furthermore, figure 4(a) in their study shows that auditors correctly detected unfulfilled requirements to report offshore wealth for 10 individuals out of the 135 individuals audited, which comes to a 7.4\% detection rate for individuals audited. Multiplying the former by the latter, gives us approximately a 0.6% detection rate. We thank Daniel Reck, one of the authors of the aforementioned study, for suggesting this approach.

# Code #
The model described is solved using a set of computer programs written in C. The code, along with the Python scripts used to process the program's output are contained in this repository. The source code is contained in the folder <a href="c/src">c/src</a>. The binary executables, which are created by compiling the program, are contained in the folder <a href="c/bin">c/bin</a>. The output of the program, which is created by running the executables, is contained in the folder <a href="c/output">c/output</a>. The makefiles used to create the executables, are contained in the top level of the <a href="c">c</a> folder. The scripts used to process the output are contained in the <a href="python">python</a> folder. The tables and figures shown in the paper and appendix are in the <a href="python/output">python/output</a> folder.

## Programs and system requirements ##
There are three main programs: `model,` `model_2type,` and `optpol`. All three programs write output in the form of CSV files. There is one file per equilibrium. The first line contains variable names (e.g. Y for GDP, WtaxRev_lost for wealth tax revenues lost to evasion...). Each row contains the values for a given period. A stationary equilibrium output file has one line. A transition has many lines, one per period.  
<br/>
**model.** This program solves for the model's equilibrium (both in the long-run and transition dynamics) for a given set of parameters. It uses OpenMP to parallelize the solution of the household's problem and iteration of the distribution. It can be run by simply typing ./bin/model from the command line. This program has many command-line options that allow the user to run the baseline and no-evasion counterfactuals, various sensitivity analyses, and transition dynamics. To see all the options, run `./bin/model --help`. Note that this program can be run in principle on any computer, but it requires a large amount of memory and lots of CPU cores to run in practice. We used a dual-CPU Xeon workstation with 40 cores and 92GB of RAM. It takes several hours to solve for a single equilibrium and at least a week to solve for a transition.

**model_2type.** This program solves the model in which some households can evade while others cannot. This version requires a different approach to memory management and the easiest way to do it was to simply create an entirely separate version of the code.

  **optpol.** This program conducts a global search for the optimal progressive wealth tax using the differential evolution algorithm (https://en.wikipedia.org/wiki/Differential_evolution). This algorithm is implemented using MPI to parallelize the solution of many steady states (associated with different tax parameters) simultaneously. This program still uses OpenMP to parallelize the household problem and distribution updating. In other words, it is a hybrid OpenMP-MPI approach. This program must be run on a supercomputer cluster such as the University of Toronto's Niagara system. We used 100 compute nodes (each with 40 cores) to run this program, assigning 4 MPI tasks to each node. Thus, each iteration of the optimization algorithm solves for 400 equilibria simultaneously. The bash script optimize.sh in the main <a href="c">c</a> folder contains the batch processing submission request we used.
  
To compile the programs, simply navigate to the <a href="c">c</a> folder and type `make model` or `make optpol` in the command line. In addition to the standard C codebase, these programs require the GNU gcc compiler (https://gcc.gnu.org) or the Intel icc compiler (https://www.intel.com/content/www/us/en/develop/documentation/cpp-compiler-developer-guide-and-reference/top.html) the GNU GSL library (https://www.gnu.org/software/gsl), OpenMP (https://www.openmp.org), and OpenMPI (https://www.open-mpi.org). We used Ubuntu Linux 20.04 to compile and run the programs, and we cannot guarantee that these programs will work without modifications in Windows or other operating systems.

## Source code ##
The source code is broken down into several files.

**calibration.c, calibration.h.** This header and source file contain all of the declarations and assignments of the model's parameters, along with several small utility routines.

**eqm.c, eqm.h.** This header and source file contain the code to solve the household's problem, update the distribution, solve for a long-run equilibrium, and solve the transition dynamics given a set of parameter values.

**main.c.** This file contains the `main(int argc, char * argv)` function for the program `model`

**main_2type.c**, **eqm_2type.c**, *eqm_2type.h.** Analogous source files for the version of the model in which some households can evade while others cannot.

**diff_evo_mpi.c, diff_evo_utils.c, diff_evo_utils.h, externs.c.** These files contain the code the for MPI implementation of the differential evolution global optimization algorithm in the program `optpol`

## Processing scripts ##
In addition to the program's source code, we use several Python scripts to create the tables and figures shown in the paper and appendix.

**tables_figs_main.py.** This script creates table 4 and figures 1-2, which show the main long-run results.

**table_dist.py.** This script creates table 3, which shows how welfare consequences of capital taxes are distributed.

**tables_figs_trans.py** This script creates figure 6 and table 7, which show transition dynamics.

**table_figs_sens.py** This script creates table 6 and figure D.1, which show the results of the sensitivity analyses we discuss in the main text, as well as versions of the model where detection revenues do not enter the government's budget constraint and hidden wealth cannot be collateralized.

**figs_chi0.py.** This script creates figures D.2-D.3, which show the effects of capital income and wealth taxes when hidden wealth cannot be collateralized

**figs_taul.py.** This script creates figures D.4-D.5, which show the effects of capital income and wealth taxes when additional revenues are used to reduce labor income taxes.

**figs_soe.py.** This script creates figures D.6-D.7, which show the results in our small-open-economy sensitivity analysis.

**figs_noconst.py.** This script creates figures D.8-D.9, which show the results in our sensitivity analysis in which there is no external financing constraint.
