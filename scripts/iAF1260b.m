fileName = 'iAF1260b.mat';
if ~exist('modelOri','var')
    modelOri = readCbModel(fileName);
end

modelOri = changeRxnBounds(modelOri,'ATPM',1000,'u');
model = modelOri;

%model.rxns

model = changeObjective(model,'BIOMASS_Ec_iAF1260_core_59p81M');
FBAsolution = optimizeCbModel(model,'max');
