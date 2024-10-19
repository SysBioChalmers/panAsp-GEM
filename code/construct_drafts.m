% replace original gene IDs with orthology cluster labels to construct
% Aspergillus_draft for each species

% Load orthology info from BPGA
BPGA=readtable('../data/genome/BPGA2ortho_GEM_custom.csv', 'Delimiter', ';');

%% Fumigatus
model=importModel("../model/templates/pan_afm.xml",false,false,true);
translationTable=id2clust(BPGA,'Afumigatus');

[draftModel,removedRxns]=getModelFromOrthology(model,translationTable);
draftModel.id='Aspergillus_fumigatus';
draftModel.name='Genome-scale model for Aspergillus fumigatus';
exportModel(draftModel,"../model/Afu_draft.xml");
exportToExcelFormat(draftModel,"../model/Afu_draft.xlsx");
%%

%% Oryzae
model=importModel("../model/templates/iWV1314.xml",false,false,true);
%%getModelFromOrthology fails w/o annotation: add placeholder
model.annotation=struct("taxonomy",'');
translationTable=id2clust(BPGA,'Aoryzae');

[draftModel,removedRxns]=getModelFromOrthology(model,translationTable);
draftModel.id='Aspergillus_oryzae';
exportModel(draftModel,"../model/Aor_draft.xml");
exportToExcelFormat(draftModel,"../model/Aor_draft.xlsx");
%%

%% Niger
model=importModel("../model/templates/iJB1325.xml",false,false,true);
%%getModelFromOrthology fails w/o annotation: add placeholder
model.annotation=struct("taxonomy",'');
translationTable=id2clust(BPGA,'Aniger');

[draftModel,removedRxns]=getModelFromOrthology(model,translationTable);
draftModel.id='Aspergillus_niger';
draftModel.name='Genome-scale metabolic model for Aspergillus niger';
exportModel(draftModel,"../model/Ani_draft.xml");
exportToExcelFormat(draftModel,"../model/Ani_draft.xlsx");
%%