clear all 

fileName = 'iMM904.mat';
if ~exist('modelOri','var')
    modelOri = readCbModel(fileName);
end

modelOri = changeRxnBounds(modelOri,'ATPM',1000,'u');
model = modelOri;

%model.rxns

model = changeObjective(model,'BIOMASS_SC5_notrace');
FBAsolution = optimizeCbModel(model,'max');
