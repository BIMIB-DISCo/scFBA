**scFBA (single-cell Flux Balance Analysis)**
==============================================================


**scFBA** is a suite of MATLAB functions to perform a single-cell Flux Balance Analysis (scFBA).


scFBA includes all the functions that are necessary to predict single-cell fluxes from the following inputs:
1) a Matlab matrix *T* including the scRNA-seq of *Ncells* single-cells, and (if available) of the bulk
2) the *SBML* of a template metabolic network *A* 


The final OUTPUT is a matrix including the flux of each reaction in each cell.


For information about the usage of each intermediate function, see the function Help.


An example of the pipeline usage is given in *scFBA_test.m*


The script *replicateScFBAresults.m* replicates all the results presented in Damiani et al, 2018.


**REFERENCE**:
 - Integration of single-cell RNA-seq data into metabolic models to characterize tumour cell populations - 
Authors: Chiara Damiani, Davide Maspero, Marzia Di Filippo, Riccardo Colombo, Dario Pescini, Alex Graudenzi, 
Hans Victor Westerhoff, Lilia Alberghina, Marco Vanoni and Giancarlo Mauri. 2018 (preprint: https://www.biorxiv.org/content/early/2018/01/30/256644)

*SYSTEM REQUIREMENTS*:
- COBRA Toolbox 2017 is required
- Gurobi solver is recommended, but scFBA works with any solver