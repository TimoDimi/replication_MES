
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Replication Material for the Paper “Regressions under Adverse Conditions”

<!-- badges: start -->
<!-- badges: end -->

The code in this repository generates the results of the simulations and
the empirical applications of the working paper Dimitriadis, T. and
Hoga, Y. (2024), Regressions under Adverse Conditions, available on
[arXiv](https://arxiv.org/abs/2311.13327).

## Data Availability

We unfortunately do not have the licence to make the bubble data from
the systemic risk application in Section 5.1 publicly available. Hence,
the code of this applications does of course not work without these
files. We still publish the code to make it available for inspection in
this repository. The data for the other two applications is included in
the repository such that these can be replicated. The simulations of
course work without any additional files.

## Simulations

The files for the simulations are available in the folder ‘simulations’.
The file ‘sim_MES.R’ generates the parameter estimates for the
regression models in Section 4, which are saved under
‘simulations/data’. It is recommended to run the file on a cluster on
10+ kernels. The file ‘sim_MES_eval.R’ evaluates these results and
generates the output in the ‘simulations/output’ folder.

## Application 1: Systemic Risk Regressions

The regression coefficients are estimated and printed using the file
‘appl1_SystemicRisk.R’. As noted above, this application cannot be
replicated with the given code as the necessary data file
‘data_bubbles.rds’ cannot be made publicly available. The respective
output used in the paper is saved to the folder ‘applications/output’

## Application 2: Managing and Allocating Portfolio Risks

The MES forecasts and (forecasted) portfolio weights are calculated in
the file ‘appl2_ERC_Portfolio.R’. The performance is then evaluated
using the file ‘appl2_ERC_Portfolio_eval.R’. All necessary data files
for this application are publicly available in this repository. The
respective output used in the paper is saved to the folder
‘applications/output’

## Application 3: Dissecting GDP Growth Vulnerabilities

The model estimates, respective output and plots are generated in the
file ‘appl3_GDP_risk.R’. All necessary data files for this application
are publicly available in this repository. The respective output used in
the paper is saved to the folder ‘applications/output’
