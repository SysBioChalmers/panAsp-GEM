% replace original gene IDs with orthology cluster labels to construct
% Aspergillus_draft for each species

% Load orthology info from BPGA
BPGA=readtable('../data/genome/BPGA2ortho_GEM_custom.csv', 'Delimiter', ';');

%% Fumigatus
model=importModel("../model/templates/pan_afm.xml",false,false,true);
translationTable=id2clust(BPGA,'Afumigatus');

[draftModel,removedRxns]=getModelFromOrthology_local(model,translationTable);

% Retain link between cluster ID and original ID in ShortNames
[~,position] = ismember(draftModel.genes,table2cell(BPGA(:,'cluster')));
draftModel.geneShortNames = table2cell(BPGA(position,'Afumigatus'));

% Remove geneMiriams field, do not link multiple miriams with same cluster
draftModel = rmfield(draftModel,"geneMiriams")

% Prepare for export
draftModel.id='Aspergillus_fumigatus';
draftModel.name='Genome-scale model for Aspergillus fumigatus';

% Export
exportModel(draftModel,"../model/Afu_draft.xml");
exportToExcelFormat(draftModel,"../model/Afu_draft.xlsx");
%%

%% Oryzae
model=importModel("../model/templates/iWV1314.xml",false,false,true);
%%getModelFromOrthology fails w/o annotation: add placeholder
model.annotation=struct("taxonomy",'');
translationTable=id2clust(BPGA,'Aoryzae');

[draftModel,removedRxns]=getModelFromOrthology_local(model,translationTable);

% For the 45 IDs starting with ZY, remove all ZYs that are isoenzymes or
% or subunits from GEM. For the remain 8 ZYs, take out GPR rule for 8 ZYs.
% For the 1 ID starting with Chr4, remove from GEM.
target=translationTable(find(contains(translationTable(:,1),("ZY"|"Chr4"))),:);
keep={'ZY006430','ZY080692','ZY088142','ZY097848', ...
    'ZY104986','ZY110993','ZY111007','ZY140875'};
genesRm=target(find(~contains(target(:,1),keep)),2);
[draftModel, affectedRxns, originalGPRs, deletedReactions]=removeGenesFromModel(draftModel, genesRm);

genesKeep=target(find(contains(target(:,1),keep)),2);
[~,position]=ismember(genesKeep,draftModel.grRules);
draftModel.grRules(position)={''};

% Retain link between cluster ID and original ID in ShortNames
[~,position] = ismember(draftModel.genes,table2cell(BPGA(:,'cluster')));
draftModel.geneShortNames = table2cell(BPGA(position,'Aoryzae'));

% Prepare for export
draftModel.id='Aspergillus_oryzae';

% Export
exportModel(draftModel,"../model/Aor_draft.xml");
exportToExcelFormat(draftModel,"../model/Aor_draft.xlsx");

%% Niger
model=importModel("../model/templates/iJB1325.xml",false,false,true);

% getModelFromOrthology fails w/o annotation: add placeholder
model.annotation=struct("taxonomy",'');
translationTable=id2clust(BPGA,'Aniger');
[draftModel,removedRxns]=getModelFromOrthology_local(model,translationTable);

% Take out GPR rule for Unknown5
genesKeep=translationTable(find(contains(translationTable(:,1),'Unknown5')),2);
[~,position]=ismember(genesKeep,draftModel.grRules);
draftModel.grRules(position)={''};

% Retain link between cluster ID and original ID in ShortNames
[~,position] = ismember(draftModel.genes,table2cell(BPGA(:,'cluster')));
draftModel.geneShortNames = table2cell(BPGA(position,'Aniger'));

% Prepare for export
draftModel.id='Aspergillus_niger';
draftModel.name='Genome-scale metabolic model for Aspergillus niger';

% Prepare for export
exportModel(draftModel,"../model/Ani_draft.xml");
exportToExcelFormat(draftModel,"../model/Ani_draft.xlsx");
%%