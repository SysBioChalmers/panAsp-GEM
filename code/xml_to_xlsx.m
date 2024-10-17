%%%%%%%%%%%%%% Aspergillus Oryzae %%%%%%%%%%%%%%

% Inspecting the Aspergillus oryzae GEM from Vongsangnak et al. (2008)
% xml obtained from:
% https://github.com/opencobra/m_model_collection/blob/master/sbml3/iWV1314.xml

model=importModel("../model/templates/iWV1314.xml",false,false,true)
% obtained without RAVEN/libSBML errors

exportToExcelFormat(model,"../model/templates/iWV1314.xlsx")

%%%%%%%%%%%%%% Aspergillus fumigatus %%%%%%%%%%%%%%
% Inspecting the Aspergillus fumigatus GEM from Mirhakkak et al. (2023)
% xml obtained from:
% https://www.ebi.ac.uk/biomodels/MODEL2211100001#Files (the pan GEM)

model=importModel("../model/templates/pan_afm.xml",false,false,true)
% obtained without RAVEN/libSBML errors

exportToExcelFormat(model,"../model/templates/pan_afm.xlsx")

%%%%%%%%%%%%%% Aspergillus niger %%%%%%%%%%%%%%
% Inspecting the Aspergillus niger GEM from Brandl et al. (2018)
% xml obtained from:
% https://fungalbiolbiotech.biomedcentral.com/articles/10.1186/s40694-018-0060-7

%%%%%%%%%%%% import
model=importModel("../model/templates/iJB1325.xml",false,false,true)

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
exportToExcelFormat(model,"../model/templates/iJB1325.xlsx")