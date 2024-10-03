% Inspecting the Aspergillus oryzae GEM from Vongsangnak et al. (2008)
% https://github.com/opencobra/m_model_collection/blob/master/sbml3/iWV1314.xml

model=importModel("iWV1314.xml",false,false,true)

% obtained without RAVEN/libSBML errors

exportToExcelFormat(model,"iWV1314.xlsx")