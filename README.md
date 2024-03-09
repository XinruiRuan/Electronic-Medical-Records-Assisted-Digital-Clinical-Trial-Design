# Electronic-Medical-Records-Assisted-Digital-Clinical-Trial-Design


This repository contains R code for conducting simulation studies and case studies. The code files are organized as follows:

1. `Simulation.R`:
   - Description: Run this code to conduct Monte Carlo simulations for both transportability and imperfect transportability scenarios.
   - Usage: Run this code to generate simulation results.

2. `plot.R`:
   - Description: Use this code to visualize the simulation results, creating Figures 1 and 2 as presented in the paper.
   - Usage: Run this code to generate plots based on simulation results.

3. `CaseStudy.R`:
   - Description: This code is used for conducting a case study, estimating treatment effects and confidence intervals.
   - Usage: Run this code to perform the HIV cash transfer case study and generate results for Figure 3.

4. `Estimators.R`:
   - Description: Contains code for data integration estimators.
   - Usage: These functions can be called from other codes for data integration purposes.

5. `Adaptive_allocation.R`:
   - Description: Includes code for solving optimal design strategies and designing adaptive experiments.
   - Usage: These functions can be used to design and implement digital clinical trials.

6. `Percentile_bootstrap.R`:
   - Description: Contains code for constructing confidence intervals using percentile bootstrap.
   - Usage: These functions can be used to construct confidence intervals.

Getting Started:
    Run the R codes in the following order: `Simulation.R`, `plot.R`, `CaseStudy.R`.
