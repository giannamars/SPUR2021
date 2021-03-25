clear all 

fileName = 'iHN637.mat';
if ~exist('modelOri','var')
    modelOri = readCbModel(fileName);
end

modelOri = changeRxnBounds(modelOri,'ATPM',1000,'u');
model = modelOri;

%model.rxns

model = changeObjective(model,'BIOMASS_Cl_DSM_WT_46p666M1');
FBAsolution = optimizeCbModel(model,'max');
