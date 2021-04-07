clear all 

fileName = 'iLJ478.mat';
if ~exist('modelOri','var')
    modelOri = readCbModel(fileName);
end

modelOri = changeRxnBounds(modelOri,'ATPM',1000,'u');
model = modelOri;

%model.rxns

model = changeObjective(model,'BIOMASS_Ecoli_TM');
FBAsolution = optimizeCbModel(model,'max');
