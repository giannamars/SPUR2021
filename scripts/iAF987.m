fileName = 'iAF987.mat';
if ~exist('modelOri','var')
    modelOri = readCbModel(fileName);
end

modelOri = changeRxnBounds(modelOri,'ATPM',1000,'u');
model = modelOri;

%model.rxns

model = changeObjective(model,'BIOMASS_Gm_GS15_core_79p20M');
FBAsolution = optimizeCbModel(model,'max');
