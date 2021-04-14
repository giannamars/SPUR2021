fileName = 'iSBO_1134.mat';
if ~exist('modelOri','var')
    modelOri = readCbModel(fileName);
end

modelOri = changeRxnBounds(modelOri,'ATPM',1000,'u');
model = modelOri;

%model.rxns

model = changeObjective(model,'BIOMASS_Ec_iJO1366_core_53p95M');
FBAsolution = optimizeCbModel(model,'max');
