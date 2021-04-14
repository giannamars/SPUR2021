fileName = 'iEK1008.mat';
if ~exist('modelOri','var')
    modelOri = readCbModel(fileName);
end

modelOri = changeRxnBounds(modelOri,'ATPM',1000,'u');
model = modelOri;

%model.rxns

model = changeObjective(model,'BIOMASS__2');
FBAsolution = optimizeCbModel(model,'max');

