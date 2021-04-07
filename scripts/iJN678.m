clear all 

fileName = 'iJN678.mat';
if ~exist('modelOri','var')
    modelOri = readCbModel(fileName);
end

modelOri = changeRxnBounds(modelOri,'ATPM',1000,'u');
model = modelOri;

%model.rxns

model = changeObjective(model,'BIOMASS_Ec_SynMixo');
FBAsolution = optimizeCbModel(model,'max');
