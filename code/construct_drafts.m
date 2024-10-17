% replace original gene IDs with orthology cluster labels

% Fumigatus
model=importModel("../model/templates/pan_afm.xml",false,false,true)
id2clust=readtable("../data/genome/fumi_temp.csv")
id2clust = table2cell(id2clust)

[draftModel,removedRxns]=getModelFromOrthology(model,id2clust)
exportModel(draftModel,"../model/Afu_draft.xml")
exportToExcelFormat(draftModel,"../model/Afu_draft.xlsx")

% Oryzae
model=importModel("../model/templates/iWV1314.xml",false,false,true)
%%getModelFromOrthology fails w/o annotation: add placeholder
model.annotation=struct("placeholder","placeholder") 
id2clust=readtable("../data/genome/ory_temp.csv")
id2clust=table2cell(id2clust)

[draftModel,removedRxns]=getModelFromOrthology(model,id2clust)

exportModel(draftModel,"../model/Aory_draft.xml")
exportToExcelFormat(draftModel,"../model/Aory_draft.xlsx")

% Niger
model=importModel("../model/templates/iJB1325.xml",false,false,true)
%%getModelFromOrthology fails w/o annotation: add placeholder
model.annotation=struct("placeholder","placeholder")

opts = detectImportOptions("../data/genome/nig_temp.csv");
opts = setvartype(opts,"Aniger",'char');  %or 'char' if you prefer
id2clust=readtable("../data/genome/nig_temp.csv",opts)
id2clust=table2cell(id2clust)

[draftModel,removedRxns]=getModelFromOrthology(model,id2clust)

exportModel(draftModel,"../model/Anig_draft.xml")
exportToExcelFormat(draftModel,"../model/Anig_draft.xlsx")