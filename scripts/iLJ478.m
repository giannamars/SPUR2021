

fileName = 'iLJ478.mat';
if ~exist('modelOri','var')
    modelOri = readCbModel(fileName);
end

modelOri = changeRxnBounds(modelOri,'ATPM',20,'u');
model = modelOri;

%model.rxns

model = changeObjective(model,'BIOMASS_Ecoli_TM');
FBAsolution = optimizeCbModel(model,'max');

%dglucose
model = changeRxnBounds(model, 'EX_glc__D_e', -18.5, 'l');

%fumarate
model = changeRxnBounds(model, 'EX_fum_e', -18.5, 'l');

%acetate
model = changeRxnBounds(model, 'EX_ac_e', -18.5, 'l');

%acetaldehyde
model = changeRxnBounds(model, 'EX_acald_e', -18.5, 'l');

%2-oxoglutarate
model = changeRxnBounds(model, 'EX_akg_e', -18.5, 'l');

%ethanol
model = changeRxnBounds(model, 'EX_etoh_e', -18.5, 'l');

%formate
model = changeRxnBounds(model, 'EX_for_e', -18.5, 'l');

%d-fructose
model = changeRxnBounds(model, 'EX_fru_e', -18.5, 'l');

%L-glutamine
model = changeRxnBounds(model, 'EX_gln__L_e', -18.5, 'l');

%L-glutamate
model = changeRxnBounds(model, 'EX_glu__L_e', -18.5, 'l');

%D-lactate
model = changeRxnBounds(model, 'EX_lac__D_e', -18.5, 'l');

%L-malate
model = changeRxnBounds(model, 'EX_mal__L_e', -18.5, 'l');

%Pyruvate
model = changeRxnBounds(model, 'EX_pyr_e', -18.5, 'l');

%Succinate
model = changeRxnBounds(model, 'EX_succ_e', -18.5, 'l');