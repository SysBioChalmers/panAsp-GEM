% Inspecting the Aspergillus niger GEM from Brandl et al. (2018)
% https://fungalbiolbiotech.biomedcentral.com/articles/10.1186/s40694-018-0060-7

%%%%%%%%%%%% import
% function model=importModel(fileName,removeExcMets,isSBML2COBRA,supressWarnings)
model=importModel("40694_2018_60_MOESM2_ESM.xml",false,false,true)

% The model contains 3741 errors.
% Assess using validator webtool: https://sbml.bioquant.uni-heidelberg.de/validator_servlet/index.jsp. 
% To me these errors seem not problematic, except for some data loss, info
% on pathway is lost in b.):
% a.	1324 errors are caused by the undesired presence of an attribute 
% "gem:genome" (all "Aspni7_Gene") in the GeneProduct object (info on genes). 
% b.	Analogously, 1910 errors are caused by the unwanted presence of 
% attribute "gem:subsystem" (e.g. conversions, Embden-Meyerhoff-Parnas 
% Pathway, Pentose-Phosphate Shunt, i.e., pathways) in the reaction object. 
% Note that the subsystem column in RXNS is indeed empty in the model. 
% c.	423 errors are caused by unwanted presence of attribute "gem:id” 
% (all different, e.g., gem:id="2ff24c08-24a5-4a91-a445-0632d9b119d2") in 
% the Or object. Note that 423 is the number of “or's” in the 
% GENE.ASSOCIATION column. 
% d.	Similarly, 81 errors were caused by the presence of attribute 
% "gem:id" in the AND object. Note that 81 is the number of “and's” in the 
% GENE.ASSOCIATION column.
% e.	Finally, 3 errors due to presemce of strict="true" in the SBML 
% namespace.

%%%%%%%%%%%% export
%function exportToExcelFormat(model,fileName,sortIds)
exportToExcelFormat(model,"40694_2018_60_MOESM2_ESM.xlsx")