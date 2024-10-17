% replace original gene IDs with orthology cluster labels

model=importModel("../data/BaseGEMs/Fumigatus_GEM.xml",false,false,true)
id2clust=readtable("../data/genome/fumi_temp.csv")
id2clust = table2cell(id2clust)

[draftModel,removedRxns]=getModelFromOrthology(model,id2clust)
exportModel(draftModel,"../model/Afu_draft.xml")
exportToExcelFormat(draftModel,"../model/Afu_draft.xlsx")


%fails
templateModel=model
orthologPairs=id2clust


% remove non-standard fields, if any
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
% identifies which GENE ASSOCIATION in RXNS are not empty

% Replace genes according to the mapped ortholog pairs, which should be in
% the defined format (an Nx2 cell array)
draftModel = templateModel;


[grRules,genes,rxnGeneMat] = replaceGrRules(draftModel.grRules,orthologPairs);
% fails: open up

grRules=draftModel.grRules
idMapping=orthologPairs

% get original list of genes from the grRules
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

% Update with modified gene fields
draftModel.grRules    = grRules;
draftModel.genes      = genes;
draftModel.rxnGeneMat = rxnGeneMat;


% find the index of non-empty grRules after replacing genes
postNonEmptyRuleInd = find(~cellfun(@isempty, draftModel.grRules));


% remove the rxns whose grRules become empty after replacement of orthologs
removedRxns = setdiff(preNonEmptyRuleInd, postNonEmptyRuleInd);
if any(removedRxns)
    draftModel = removeReactions(draftModel, removedRxns, true, true, true);
end