% Inspecting the Aspergillus fumigatus GEM from Mirhakkak et al. (2023)
% https://www.ebi.ac.uk/biomodels/MODEL2211100001#Files (the pan GEM)

model=importModel("Pan_Aspergillus_fumigatus.xml",false,false,true)

% obtained without RAVEN/libSBML errors

exportToExcelFormat(model,"Pan_Aspergillus_fumigatus.xlsx")

pairs = readtable("orthologPairs_Fumigatus_test.csv", "FileType","text",'Delimiter', ';');

%[draftModel, removedRxns]=getModelFromOrthology(model,pairs)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
templateModel=model
orthologPairs=pairs

fieldsToRemove = intersect({'rxnFrom','metFrom'}, fieldnames(templateModel));
templateModel = rmfield(templateModel, fieldsToRemove);

% clean metadata fields
templateModel.id = '';
templateModel.description = '';
templateModel.version = '';
templateModel.annotation = structfun(@(x) '',templateModel.annotation,'UniformOutput',0);
templateModel.annotation.defaultLB = -1000;
templateModel.annotation.defaultUB = 1000;

% find the index of non-empty grRules before replacing genes
preNonEmptyRuleInd = find(~cellfun(@isempty, templateModel.grRules));


% Replace genes according to the mapped ortholog pairs, which should be in
% the defined format (an Nx2 cell array)
draftModel = templateModel;

%[grRules,genes,rxnGeneMat] = replaceGrRules(draftModel.grRules,orthologPairs);

grRules=draftModel.grRules
idMapping=orthologPairs

genes_orig = getGenesFromGrRules(grRules);
if ismember('and',genes_orig) || ismember('or',genes_orig)
    error('Problem reading grRules. Verify that all "and" and "or" elements are lowercase and surrounded by spaces.');
end

% remove duplicate genes
genes_orig  = unique(genes_orig,'stable');

% idMapping should be a NX2 cell array
if ~(size(idMapping,2) == 2)
    error('The idMapping data structure must be a NX2 cell array.');
% input genes should be consistent with those retrieved from grRules
elseif all(~ismember(idMapping(:,1), genes_orig))
    error('The input genes are Not consistent with those retrieved from grRules.');
end

