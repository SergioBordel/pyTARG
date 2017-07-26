# pyTARG
Python library fo work with genome scale metabolic models
pyTARG is a library that contains functions to work with Genome Scale Metabolic Models with the goal of finding drug targets against cancer
1. fullconstrain(model,expr,coefficient) takes as inputs a model imported with COBRApy, a dictionary with expression levels of each gene (the keys are ENSEMBL
gene identifiers and the values are expression levels in RPMK, the coefficient takes a deffault value of 0.1. The function outputs a constrained model
2. flux(model) takes a constrained model as input and outputs a flux distribution in mmol/h/g-biomass
3. block(model,targets) takes as inputs a constrained model and a list of targets, which can be genes, reactions or metabolites. The output is the relative value of
the objective function after constraining the fluxes in the target reactions to 0.1 times their initial value
4. personal(model1,model2) takes two constrained models as inputs, model1 is the cancer model to be targetted and model2 a reference healthy cell type. The output
is a list of reactions to be targetted simultaneously.
