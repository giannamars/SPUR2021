clear all 

fileName = 'iSB619.mat';
if ~exist('modelOri','var')
    modelOri = readCbModel(fileName);
end

modelOri = changeRxnBounds(modelOri,'ATPM',1000,'u');
model = modelOri;

%model.rxns

model = changeObjective(model,'BIOMASS_SA_1a');
FBAsolution = optimizeCbModel(model,'max');
