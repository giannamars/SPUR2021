fileName = 'iEK1008.mat';
if ~exist('modelOri','var')
    modelOri = readCbModel(fileName);
end

modelOri = changeRxnBounds(modelOri,'ATPM',1000,'u');
model = modelOri;

%model.rxns
%model = changeObjective(model,'BIOMASS__2');
%FBAsolution = optimizeCbModel(model,'max');
outmodel = writeCbModel(model, 'xls', 'iEK1008_model.xls'); % write all reactions to xls file
modelClosed = model;
% Find all exchange reactions in the model
modelexchanges1 = strmatch('Ex_', model.rxns); 
modelexchanges4 = strmatch('EX_', model.rxns);          % identify 70 exchange reactions
modelexchanges2 = strmatch('DM_', model.rxns);
modelexchanges3 = strmatch('sink_', model.rxns);
% To make sure we've found all exchange reactions via string matching, we
% also check for reactions that contain only one non-zero entry in the S matrix column
selExc = (find(full((sum(abs(model.S)==1, 1)==1) & (sum(model.S~=0) == 1))))'; % identify 70 exchange reactions
modelexchanges = unique([modelexchanges1; modelexchanges2; modelexchanges3; modelexchanges4; selExc]);
reaction_abbrvs = model.rxns(modelexchanges); % store reaction abbreviations in model
reaction_names = model.rxnNames(modelexchanges); % store reaction names

% set lower bounds (=uptake) of all exchange reactions to zero, meaning to uptake
model.lb(find(ismember(model.rxns, model.rxns(modelexchanges))))=0;
% set uber bounds (=excretion) to 1000 just to be safe
model.ub(selExc) = 1000; 
% change objective function to biomass optimization
model = changeObjective(model,'BIOMASS_');
modelClosedOri = modelClosed; % store this setup in modelClosedOri so we don't have to repeat above lines

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Acetate
modelClosed = modelClosedOri;
modelClosed.lb(find(ismember(modelClosed.rxns, 'EX_o2_e'))) = -1000;
modelClosed.lb(find(ismember(modelClosed.rxns, 'EX_h2o_e'))) = -1000;
modelClosed.ub(find(ismember(modelClosed.rxns, 'EX_h2o_e'))) = 1000;
modelClosed.lb(find(ismember(modelClosed.rxns, 'EX_co2_e'))) = -1000;
modelClosed.ub(find(ismember(modelClosed.rxns, 'EX_co2_e'))) = 1000;
modelClosed.lb(find(ismember(modelClosed.rxns, 'EX_ac_e'))) = -1;
modelClosed.ub(find(ismember(modelClosed.rxns, 'EX_ac_e'))) = -1;
FBA = optimizeCbModel(modelClosed, 'max');
% calc C
%printFluxVector(model,FBA.x,true,true) % only print exchange reaction fluxes 
%
exchange_fluxes = FBA.x(find(ismember(modelClosed.rxns, reaction_abbrvs)));
%
exchange_atoms = zeros(length(exchange_fluxes),1);
for i = 1:length(exchange_atoms)
    tmp1 = printRxnFormula(modelClosed,reaction_abbrvs(i));
    tmp2 = strtrim(regexp(tmp1{1}, '->|<=>', 'split'));
    metID = findMetIDs(modelClosed, tmp2{1});
    metFormula= modelClosed.metFormulas(metID);
    atomicStructure = parse_formula(metFormula{1});
    C_atoms = zeros(1);
    if(isfield(atomicStructure, 'C'))
        C_atoms(1) = atomicStructure.C;
    else
        C_atoms(1) = 0;
    end
    exchange_atoms(i) = C_atoms(1);
end
%
Secretion_fluxes = exchange_fluxes(exchange_fluxes > 0.0);
Secretion_atoms = exchange_atoms(exchange_fluxes > 0.0);
Uptake_fluxes = exchange_fluxes(exchange_fluxes < 0.0);
Uptake_atoms = exchange_atoms(exchange_fluxes < 0.0);
%
Uptake_tmp = abs(sum(Uptake_fluxes.*Uptake_atoms));
Secretion_tmp = abs(sum(Secretion_fluxes.*Secretion_atoms));
%
iEK1008_1_C = (Uptake_tmp - Secretion_tmp)/Uptake_tmp;   % not the right value...
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Pyruvate
modelClosed = modelClosedOri;
modelClosed.lb(find(ismember(modelClosed.rxns, 'EX_o2_e'))) = -1000;
modelClosed.lb(find(ismember(modelClosed.rxns, 'EX_h2o_e'))) = -1000;
modelClosed.ub(find(ismember(modelClosed.rxns, 'EX_h2o_e'))) = 1000;
modelClosed.lb(find(ismember(modelClosed.rxns, 'EX_co2_e'))) = -1000;
modelClosed.ub(find(ismember(modelClosed.rxns, 'EX_co2_e'))) = 1000;
modelClosed.lb(find(ismember(modelClosed.rxns, 'EX_pyr_e'))) = -1;
modelClosed.ub(find(ismember(modelClosed.rxns, 'EX_pyr_e'))) = -1;
FBA = optimizeCbModel(modelClosed, 'max');

exchange_fluxes = FBA.x(find(ismember(modelClosed.rxns, reaction_abbrvs)));%
%
exchange_atoms = zeros(length(exchange_fluxes),1);
for i = 1:length(exchange_atoms)
    tmp1 = printRxnFormula(modelClosed,reaction_abbrvs(i));
    tmp2 = strtrim(regexp(tmp1{1}, '->|<=>', 'split'));
    metID = findMetIDs(modelClosed, tmp2{1});
    metFormula= modelClosed.metFormulas(metID);
    atomicStructure = parse_formula(metFormula{1});
    C_atoms = zeros(1);
    if(isfield(atomicStructure, 'C'))
        C_atoms(1) = atomicStructure.C;
    else
        C_atoms(1) = 0;
    end
    exchange_atoms(i) = C_atoms(1);
end
%
Secretion_fluxes = exchange_fluxes(exchange_fluxes > 0.0);
Secretion_atoms = exchange_atoms(exchange_fluxes > 0.0);
Uptake_fluxes = exchange_fluxes(exchange_fluxes < 0.0);
Uptake_atoms = exchange_atoms(exchange_fluxes < 0.0);
%
Uptake_tmp = abs(sum(Uptake_fluxes.*Uptake_atoms));
Secretion_tmp = abs(sum(Secretion_fluxes.*Secretion_atoms));
%
iEK1008_2_C = (Uptake_tmp - Secretion_tmp)/Uptake_tmp;   % not the right value...
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% D-glucose
modelClosed = modelClosedOri;
modelClosed.lb(find(ismember(modelClosed.rxns, 'EX_o2_e'))) = -1000;
modelClosed.lb(find(ismember(modelClosed.rxns, 'EX_h2o_e'))) = -1000;
modelClosed.ub(find(ismember(modelClosed.rxns, 'EX_h2o_e'))) = 1000;
modelClosed.lb(find(ismember(modelClosed.rxns, 'EX_co2_e'))) = -1000;
modelClosed.ub(find(ismember(modelClosed.rxns, 'EX_co2_e'))) = 1000;
modelClosed.lb(find(ismember(modelClosed.rxns, 'EX_glc__D_e'))) = -1;
modelClosed.ub(find(ismember(modelClosed.rxns, 'EX_glc__D_e'))) = -1;
FBA = optimizeCbModel(modelClosed, 'max');

exchange_fluxes = FBA.x(find(ismember(modelClosed.rxns, reaction_abbrvs)));%
%
exchange_atoms = zeros(length(exchange_fluxes),1);
for i = 1:length(exchange_atoms)
    tmp1 = printRxnFormula(modelClosed,reaction_abbrvs(i));
    tmp2 = strtrim(regexp(tmp1{1}, '->|<=>', 'split'));
    metID = findMetIDs(modelClosed, tmp2{1});
    metFormula= modelClosed.metFormulas(metID);
    atomicStructure = parse_formula(metFormula{1});
    C_atoms = zeros(1);
    if(isfield(atomicStructure, 'C'))
        C_atoms(1) = atomicStructure.C;
    else
        C_atoms(1) = 0;
    end
    exchange_atoms(i) = C_atoms(1);
end
%
Secretion_fluxes = exchange_fluxes(exchange_fluxes > 0.0);
Secretion_atoms = exchange_atoms(exchange_fluxes > 0.0);
Uptake_fluxes = exchange_fluxes(exchange_fluxes < 0.0);
Uptake_atoms = exchange_atoms(exchange_fluxes < 0.0);
%
Uptake_tmp = abs(sum(Uptake_fluxes.*Uptake_atoms));
Secretion_tmp = abs(sum(Secretion_fluxes.*Secretion_atoms));
%
iEK1008_3_C = (Uptake_tmp - Secretion_tmp)/Uptake_tmp;   % not the right value...
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fumarate
modelClosed = modelClosedOri;
modelClosed.lb(find(ismember(modelClosed.rxns, 'EX_o2_e'))) = -1000;
modelClosed.lb(find(ismember(modelClosed.rxns, 'EX_h2o_e'))) = -1000;
modelClosed.ub(find(ismember(modelClosed.rxns, 'EX_h2o_e'))) = 1000;
modelClosed.lb(find(ismember(modelClosed.rxns, 'EX_co2_e'))) = -1000;
modelClosed.ub(find(ismember(modelClosed.rxns, 'EX_co2_e'))) = 1000;
modelClosed.lb(find(ismember(modelClosed.rxns, 'EX_fum_e'))) = -1;
modelClosed.ub(find(ismember(modelClosed.rxns, 'EX_fum_e'))) = -1;
FBA = optimizeCbModel(modelClosed, 'max');

exchange_fluxes = FBA.x(find(ismember(modelClosed.rxns, reaction_abbrvs)));%
%
exchange_atoms = zeros(length(exchange_fluxes),1);
for i = 1:length(exchange_atoms)
    tmp1 = printRxnFormula(modelClosed,reaction_abbrvs(i));
    tmp2 = strtrim(regexp(tmp1{1}, '->|<=>', 'split'));
    metID = findMetIDs(modelClosed, tmp2{1});
    metFormula= modelClosed.metFormulas(metID);
    atomicStructure = parse_formula(metFormula{1});
    C_atoms = zeros(1);
    if(isfield(atomicStructure, 'C'))
        C_atoms(1) = atomicStructure.C;
    else
        C_atoms(1) = 0;
    end
    exchange_atoms(i) = C_atoms(1);
end
%
Secretion_fluxes = exchange_fluxes(exchange_fluxes > 0.0);
Secretion_atoms = exchange_atoms(exchange_fluxes > 0.0);
Uptake_fluxes = exchange_fluxes(exchange_fluxes < 0.0);
Uptake_atoms = exchange_atoms(exchange_fluxes < 0.0);
%
Uptake_tmp = abs(sum(Uptake_fluxes.*Uptake_atoms));
Secretion_tmp = abs(sum(Secretion_fluxes.*Secretion_atoms));
%
iEK1008_4_C = (Uptake_tmp - Secretion_tmp)/Uptake_tmp;  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Acetaldehyde
modelClosed = modelClosedOri;
modelClosed.lb(find(ismember(modelClosed.rxns, 'EX_o2_e'))) = -1000;
modelClosed.lb(find(ismember(modelClosed.rxns, 'EX_h2o_e'))) = -1000;
modelClosed.ub(find(ismember(modelClosed.rxns, 'EX_h2o_e'))) = 1000;
modelClosed.lb(find(ismember(modelClosed.rxns, 'EX_co2_e'))) = -1000;
modelClosed.ub(find(ismember(modelClosed.rxns, 'EX_co2_e'))) = 1000;
modelClosed.lb(find(ismember(modelClosed.rxns, 'EX_acald_e'))) = -1;
modelClosed.ub(find(ismember(modelClosed.rxns, 'EX_acald_e'))) = -1;
FBA = optimizeCbModel(modelClosed, 'max');

exchange_fluxes = FBA.x(find(ismember(modelClosed.rxns, reaction_abbrvs)));%
%
exchange_atoms = zeros(length(exchange_fluxes),1);
for i = 1:length(exchange_atoms)
    tmp1 = printRxnFormula(modelClosed,reaction_abbrvs(i));
    tmp2 = strtrim(regexp(tmp1{1}, '->|<=>', 'split'));
    metID = findMetIDs(modelClosed, tmp2{1});
    metFormula= modelClosed.metFormulas(metID);
    atomicStructure = parse_formula(metFormula{1});
    C_atoms = zeros(1);
    if(isfield(atomicStructure, 'C'))
        C_atoms(1) = atomicStructure.C;
    else
        C_atoms(1) = 0;
    end
    exchange_atoms(i) = C_atoms(1);
end
%
Secretion_fluxes = exchange_fluxes(exchange_fluxes > 0.0);
Secretion_atoms = exchange_atoms(exchange_fluxes > 0.0);
Uptake_fluxes = exchange_fluxes(exchange_fluxes < 0.0);
Uptake_atoms = exchange_atoms(exchange_fluxes < 0.0);
%
Uptake_tmp = abs(sum(Uptake_fluxes.*Uptake_atoms));
Secretion_tmp = abs(sum(Secretion_fluxes.*Secretion_atoms));
%
iEK1008_5_C = (Uptake_tmp - Secretion_tmp)/Uptake_tmp;   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2-oxoglutarate
modelClosed = modelClosedOri;
modelClosed.lb(find(ismember(modelClosed.rxns, 'EX_o2_e'))) = -1000;
modelClosed.lb(find(ismember(modelClosed.rxns, 'EX_h2o_e'))) = -1000;
modelClosed.ub(find(ismember(modelClosed.rxns, 'EX_h2o_e'))) = 1000;
modelClosed.lb(find(ismember(modelClosed.rxns, 'EX_co2_e'))) = -1000;
modelClosed.ub(find(ismember(modelClosed.rxns, 'EX_co2_e'))) = 1000;
modelClosed.lb(find(ismember(modelClosed.rxns, 'EX_akg_e'))) = -1;
modelClosed.ub(find(ismember(modelClosed.rxns, 'EX_akg_e'))) = -1;
FBA = optimizeCbModel(modelClosed, 'max');

exchange_fluxes = FBA.x(find(ismember(modelClosed.rxns, reaction_abbrvs)));%
%
exchange_atoms = zeros(length(exchange_fluxes),1);
for i = 1:length(exchange_atoms)
    tmp1 = printRxnFormula(modelClosed,reaction_abbrvs(i));
    tmp2 = strtrim(regexp(tmp1{1}, '->|<=>', 'split'));
    metID = findMetIDs(modelClosed, tmp2{1});
    metFormula= modelClosed.metFormulas(metID);
    atomicStructure = parse_formula(metFormula{1});
    C_atoms = zeros(1);
    if(isfield(atomicStructure, 'C'))
        C_atoms(1) = atomicStructure.C;
    else
        C_atoms(1) = 0;
    end
    exchange_atoms(i) = C_atoms(1);
end
%
Secretion_fluxes = exchange_fluxes(exchange_fluxes > 0.0);
Secretion_atoms = exchange_atoms(exchange_fluxes > 0.0);
Uptake_fluxes = exchange_fluxes(exchange_fluxes < 0.0);
Uptake_atoms = exchange_atoms(exchange_fluxes < 0.0);
%
Uptake_tmp = abs(sum(Uptake_fluxes.*Uptake_atoms));
Secretion_tmp = abs(sum(Secretion_fluxes.*Secretion_atoms));
%
iEK1008_6_C = (Uptake_tmp - Secretion_tmp)/Uptake_tmp;   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ethanol
modelClosed = modelClosedOri;
modelClosed.lb(find(ismember(modelClosed.rxns, 'EX_o2_e'))) = -1000;
modelClosed.lb(find(ismember(modelClosed.rxns, 'EX_h2o_e'))) = -1000;
modelClosed.ub(find(ismember(modelClosed.rxns, 'EX_h2o_e'))) = 1000;
modelClosed.lb(find(ismember(modelClosed.rxns, 'EX_co2_e'))) = -1000;
modelClosed.ub(find(ismember(modelClosed.rxns, 'EX_co2_e'))) = 1000;
modelClosed.lb(find(ismember(modelClosed.rxns, 'EX_etoh_e'))) = -1;
modelClosed.ub(find(ismember(modelClosed.rxns, 'EX_etoh_e'))) = -1;
FBA = optimizeCbModel(modelClosed, 'max');

exchange_fluxes = FBA.x(find(ismember(modelClosed.rxns, reaction_abbrvs)));%
%
exchange_atoms = zeros(length(exchange_fluxes),1);
for i = 1:length(exchange_atoms)
    tmp1 = printRxnFormula(modelClosed,reaction_abbrvs(i));
    tmp2 = strtrim(regexp(tmp1{1}, '->|<=>', 'split'));
    metID = findMetIDs(modelClosed, tmp2{1});
    metFormula= modelClosed.metFormulas(metID);
    atomicStructure = parse_formula(metFormula{1});
    C_atoms = zeros(1);
    if(isfield(atomicStructure, 'C'))
        C_atoms(1) = atomicStructure.C;
    else
        C_atoms(1) = 0;
    end
    exchange_atoms(i) = C_atoms(1);
end
%
Secretion_fluxes = exchange_fluxes(exchange_fluxes > 0.0);
Secretion_atoms = exchange_atoms(exchange_fluxes > 0.0);
Uptake_fluxes = exchange_fluxes(exchange_fluxes < 0.0);
Uptake_atoms = exchange_atoms(exchange_fluxes < 0.0);
%
Uptake_tmp = abs(sum(Uptake_fluxes.*Uptake_atoms));
Secretion_tmp = abs(sum(Secretion_fluxes.*Secretion_atoms));
%
iEK1008_7_C = (Uptake_tmp - Secretion_tmp)/Uptake_tmp;  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Formate
modelClosed = modelClosedOri;
modelClosed.lb(find(ismember(modelClosed.rxns, 'EX_o2_e'))) = -1000;
modelClosed.lb(find(ismember(modelClosed.rxns, 'EX_h2o_e'))) = -1000;
modelClosed.ub(find(ismember(modelClosed.rxns, 'EX_h2o_e'))) = 1000;
modelClosed.lb(find(ismember(modelClosed.rxns, 'EX_co2_e'))) = -1000;
modelClosed.ub(find(ismember(modelClosed.rxns, 'EX_co2_e'))) = 1000;
modelClosed.lb(find(ismember(modelClosed.rxns, 'EX_for_e'))) = -1;
modelClosed.ub(find(ismember(modelClosed.rxns, 'EX_for_e'))) = -1;
FBA = optimizeCbModel(modelClosed, 'max');

exchange_fluxes = FBA.x(find(ismember(modelClosed.rxns, reaction_abbrvs)));%
%
exchange_atoms = zeros(length(exchange_fluxes),1);
for i = 1:length(exchange_atoms)
    tmp1 = printRxnFormula(modelClosed,reaction_abbrvs(i));
    tmp2 = strtrim(regexp(tmp1{1}, '->|<=>', 'split'));
    metID = findMetIDs(modelClosed, tmp2{1});
    metFormula= modelClosed.metFormulas(metID);
    atomicStructure = parse_formula(metFormula{1});
    C_atoms = zeros(1);
    if(isfield(atomicStructure, 'C'))
        C_atoms(1) = atomicStructure.C;
    else
        C_atoms(1) = 0;
    end
    exchange_atoms(i) = C_atoms(1);
end
%
Secretion_fluxes = exchange_fluxes(exchange_fluxes > 0.0);
Secretion_atoms = exchange_atoms(exchange_fluxes > 0.0);
Uptake_fluxes = exchange_fluxes(exchange_fluxes < 0.0);
Uptake_atoms = exchange_atoms(exchange_fluxes < 0.0);
%
Uptake_tmp = abs(sum(Uptake_fluxes.*Uptake_atoms));
Secretion_tmp = abs(sum(Secretion_fluxes.*Secretion_atoms));
%
iEK1008_8_C = (Uptake_tmp - Secretion_tmp)/Uptake_tmp;  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% D-fructose
modelClosed = modelClosedOri;
modelClosed.lb(find(ismember(modelClosed.rxns, 'EX_o2_e'))) = -1000;
modelClosed.lb(find(ismember(modelClosed.rxns, 'EX_h2o_e'))) = -1000;
modelClosed.ub(find(ismember(modelClosed.rxns, 'EX_h2o_e'))) = 1000;
modelClosed.lb(find(ismember(modelClosed.rxns, 'EX_co2_e'))) = -1000;
modelClosed.ub(find(ismember(modelClosed.rxns, 'EX_co2_e'))) = 1000;
modelClosed.lb(find(ismember(modelClosed.rxns, 'EX_fru_e'))) = -1;
modelClosed.ub(find(ismember(modelClosed.rxns, 'EX_fru_e'))) = -1;
FBA = optimizeCbModel(modelClosed, 'max');

exchange_fluxes = FBA.x(find(ismember(modelClosed.rxns, reaction_abbrvs)));%
%
exchange_atoms = zeros(length(exchange_fluxes),1);
for i = 1:length(exchange_atoms)
    tmp1 = printRxnFormula(modelClosed,reaction_abbrvs(i));
    tmp2 = strtrim(regexp(tmp1{1}, '->|<=>', 'split'));
    metID = findMetIDs(modelClosed, tmp2{1});
    metFormula= modelClosed.metFormulas(metID);
    atomicStructure = parse_formula(metFormula{1});
    C_atoms = zeros(1);
    if(isfield(atomicStructure, 'C'))
        C_atoms(1) = atomicStructure.C;
    else
        C_atoms(1) = 0;
    end
    exchange_atoms(i) = C_atoms(1);
end
%
Secretion_fluxes = exchange_fluxes(exchange_fluxes > 0.0);
Secretion_atoms = exchange_atoms(exchange_fluxes > 0.0);
Uptake_fluxes = exchange_fluxes(exchange_fluxes < 0.0);
Uptake_atoms = exchange_atoms(exchange_fluxes < 0.0);
%
Uptake_tmp = abs(sum(Uptake_fluxes.*Uptake_atoms));
Secretion_tmp = abs(sum(Secretion_fluxes.*Secretion_atoms));
%
iEK1008_9_C = (Uptake_tmp - Secretion_tmp)/Uptake_tmp;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% L-glutamine
modelClosed = modelClosedOri;
modelClosed.lb(find(ismember(modelClosed.rxns, 'EX_o2_e'))) = -1000;
modelClosed.lb(find(ismember(modelClosed.rxns, 'EX_h2o_e'))) = -1000;
modelClosed.ub(find(ismember(modelClosed.rxns, 'EX_h2o_e'))) = 1000;
modelClosed.lb(find(ismember(modelClosed.rxns, 'EX_co2_e'))) = -1000;
modelClosed.ub(find(ismember(modelClosed.rxns, 'EX_co2_e'))) = 1000;
modelClosed.lb(find(ismember(modelClosed.rxns, 'EX_gln__L_e'))) = -1;
modelClosed.ub(find(ismember(modelClosed.rxns, 'EX_gln__L_e'))) = -1;
FBA = optimizeCbModel(modelClosed, 'max');

exchange_fluxes = FBA.x(find(ismember(modelClosed.rxns, reaction_abbrvs)));%
%
exchange_atoms = zeros(length(exchange_fluxes),1);
for i = 1:length(exchange_atoms)
    tmp1 = printRxnFormula(modelClosed,reaction_abbrvs(i));
    tmp2 = strtrim(regexp(tmp1{1}, '->|<=>', 'split'));
    metID = findMetIDs(modelClosed, tmp2{1});
    metFormula= modelClosed.metFormulas(metID);
    atomicStructure = parse_formula(metFormula{1});
    C_atoms = zeros(1);
    if(isfield(atomicStructure, 'C'))
        C_atoms(1) = atomicStructure.C;
    else
        C_atoms(1) = 0;
    end
    exchange_atoms(i) = C_atoms(1);
end
%
Secretion_fluxes = exchange_fluxes(exchange_fluxes > 0.0);
Secretion_atoms = exchange_atoms(exchange_fluxes > 0.0);
Uptake_fluxes = exchange_fluxes(exchange_fluxes < 0.0);
Uptake_atoms = exchange_atoms(exchange_fluxes < 0.0);
%
Uptake_tmp = abs(sum(Uptake_fluxes.*Uptake_atoms));
Secretion_tmp = abs(sum(Secretion_fluxes.*Secretion_atoms));
%
iEK1008_10_C = (Uptake_tmp - Secretion_tmp)/Uptake_tmp;  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% L-glutamate
modelClosed = modelClosedOri;
modelClosed.lb(find(ismember(modelClosed.rxns, 'EX_o2_e'))) = -1000;
modelClosed.lb(find(ismember(modelClosed.rxns, 'EX_h2o_e'))) = -1000;
modelClosed.ub(find(ismember(modelClosed.rxns, 'EX_h2o_e'))) = 1000;
modelClosed.lb(find(ismember(modelClosed.rxns, 'EX_co2_e'))) = -1000;
modelClosed.ub(find(ismember(modelClosed.rxns, 'EX_co2_e'))) = 1000;
modelClosed.lb(find(ismember(modelClosed.rxns, 'EX_glu__L_e'))) = -1;
modelClosed.ub(find(ismember(modelClosed.rxns, 'EX_glu__L_e'))) = -1;
FBA = optimizeCbModel(modelClosed, 'max');

exchange_fluxes = FBA.x(find(ismember(modelClosed.rxns, reaction_abbrvs)));%
%
exchange_atoms = zeros(length(exchange_fluxes),1);
for i = 1:length(exchange_atoms)
    tmp1 = printRxnFormula(modelClosed,reaction_abbrvs(i));
    tmp2 = strtrim(regexp(tmp1{1}, '->|<=>', 'split'));
    metID = findMetIDs(modelClosed, tmp2{1});
    metFormula= modelClosed.metFormulas(metID);
    atomicStructure = parse_formula(metFormula{1});
    C_atoms = zeros(1);
    if(isfield(atomicStructure, 'C'))
        C_atoms(1) = atomicStructure.C;
    else
        C_atoms(1) = 0;
    end
    exchange_atoms(i) = C_atoms(1);
end
%
Secretion_fluxes = exchange_fluxes(exchange_fluxes > 0.0);
Secretion_atoms = exchange_atoms(exchange_fluxes > 0.0);
Uptake_fluxes = exchange_fluxes(exchange_fluxes < 0.0);
Uptake_atoms = exchange_atoms(exchange_fluxes < 0.0);
%
Uptake_tmp = abs(sum(Uptake_fluxes.*Uptake_atoms));
Secretion_tmp = abs(sum(Secretion_fluxes.*Secretion_atoms));
%
iEK1008_11_C = (Uptake_tmp - Secretion_tmp)/Uptake_tmp; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% D-lactate
modelClosed = modelClosedOri;
modelClosed.lb(find(ismember(modelClosed.rxns, 'EX_o2_e'))) = -1000;
modelClosed.lb(find(ismember(modelClosed.rxns, 'EX_h2o_e'))) = -1000;
modelClosed.ub(find(ismember(modelClosed.rxns, 'EX_h2o_e'))) = 1000;
modelClosed.lb(find(ismember(modelClosed.rxns, 'EX_co2_e'))) = -1000;
modelClosed.ub(find(ismember(modelClosed.rxns, 'EX_co2_e'))) = 1000;
modelClosed.lb(find(ismember(modelClosed.rxns, 'EX_lac__D_e'))) = -1;
modelClosed.ub(find(ismember(modelClosed.rxns, 'EX_lac__D_e'))) = -1;
FBA = optimizeCbModel(modelClosed, 'max');

exchange_fluxes = FBA.x(find(ismember(modelClosed.rxns, reaction_abbrvs)));%
%
exchange_atoms = zeros(length(exchange_fluxes),1);
for i = 1:length(exchange_atoms)
    tmp1 = printRxnFormula(modelClosed,reaction_abbrvs(i));
    tmp2 = strtrim(regexp(tmp1{1}, '->|<=>', 'split'));
    metID = findMetIDs(modelClosed, tmp2{1});
    metFormula= modelClosed.metFormulas(metID);
    atomicStructure = parse_formula(metFormula{1});
    C_atoms = zeros(1);
    if(isfield(atomicStructure, 'C'))
        C_atoms(1) = atomicStructure.C;
    else
        C_atoms(1) = 0;
    end
    exchange_atoms(i) = C_atoms(1);
end
%
Secretion_fluxes = exchange_fluxes(exchange_fluxes > 0.0);
Secretion_atoms = exchange_atoms(exchange_fluxes > 0.0);
Uptake_fluxes = exchange_fluxes(exchange_fluxes < 0.0);
Uptake_atoms = exchange_atoms(exchange_fluxes < 0.0);
%
Uptake_tmp = abs(sum(Uptake_fluxes.*Uptake_atoms));
Secretion_tmp = abs(sum(Secretion_fluxes.*Secretion_atoms));
%
iEK1008_12_C = (Uptake_tmp - Secretion_tmp)/Uptake_tmp; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% L-malate
modelClosed = modelClosedOri;
modelClosed.lb(find(ismember(modelClosed.rxns, 'EX_o2_e'))) = -1000;
modelClosed.lb(find(ismember(modelClosed.rxns, 'EX_h2o_e'))) = -1000;
modelClosed.ub(find(ismember(modelClosed.rxns, 'EX_h2o_e'))) = 1000;
modelClosed.lb(find(ismember(modelClosed.rxns, 'EX_co2_e'))) = -1000;
modelClosed.ub(find(ismember(modelClosed.rxns, 'EX_co2_e'))) = 1000;
modelClosed.lb(find(ismember(modelClosed.rxns, 'EX_mal__L_e'))) = -1;
modelClosed.ub(find(ismember(modelClosed.rxns, 'EX_mal__L_e'))) = -1;
FBA = optimizeCbModel(modelClosed, 'max');

exchange_fluxes = FBA.x(find(ismember(modelClosed.rxns, reaction_abbrvs)));%
%
exchange_atoms = zeros(length(exchange_fluxes),1);
for i = 1:length(exchange_atoms)
    tmp1 = printRxnFormula(modelClosed,reaction_abbrvs(i));
    tmp2 = strtrim(regexp(tmp1{1}, '->|<=>', 'split'));
    metID = findMetIDs(modelClosed, tmp2{1});
    metFormula= modelClosed.metFormulas(metID);
    atomicStructure = parse_formula(metFormula{1});
    C_atoms = zeros(1);
    if(isfield(atomicStructure, 'C'))
        C_atoms(1) = atomicStructure.C;
    else
        C_atoms(1) = 0;
    end
    exchange_atoms(i) = C_atoms(1);
end
%
Secretion_fluxes = exchange_fluxes(exchange_fluxes > 0.0);
Secretion_atoms = exchange_atoms(exchange_fluxes > 0.0);
Uptake_fluxes = exchange_fluxes(exchange_fluxes < 0.0);
Uptake_atoms = exchange_atoms(exchange_fluxes < 0.0);
%
Uptake_tmp = abs(sum(Uptake_fluxes.*Uptake_atoms));
Secretion_tmp = abs(sum(Secretion_fluxes.*Secretion_atoms));
%
iEK1008_13_C = (Uptake_tmp - Secretion_tmp)/Uptake_tmp; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Succinate
modelClosed = modelClosedOri;
modelClosed.lb(find(ismember(modelClosed.rxns, 'EX_o2_e'))) = -1000;
modelClosed.lb(find(ismember(modelClosed.rxns, 'EX_h2o_e'))) = -1000;
modelClosed.ub(find(ismember(modelClosed.rxns, 'EX_h2o_e'))) = 1000;
modelClosed.lb(find(ismember(modelClosed.rxns, 'EX_co2_e'))) = -1000;
modelClosed.ub(find(ismember(modelClosed.rxns, 'EX_co2_e'))) = 1000;
modelClosed.lb(find(ismember(modelClosed.rxns, 'EX_succ_e'))) = -1;
modelClosed.ub(find(ismember(modelClosed.rxns, 'EX_succ_e'))) = -1;
FBA = optimizeCbModel(modelClosed, 'max');

exchange_fluxes = FBA.x(find(ismember(modelClosed.rxns, reaction_abbrvs)));%
%
exchange_atoms = zeros(length(exchange_fluxes),1);
for i = 1:length(exchange_atoms)
    tmp1 = printRxnFormula(modelClosed,reaction_abbrvs(i));
    tmp2 = strtrim(regexp(tmp1{1}, '->|<=>', 'split'));
    metID = findMetIDs(modelClosed, tmp2{1});
    metFormula= modelClosed.metFormulas(metID);
    atomicStructure = parse_formula(metFormula{1});
    C_atoms = zeros(1);
    if(isfield(atomicStructure, 'C'))
        C_atoms(1) = atomicStructure.C;
    else
        C_atoms(1) = 0;
    end
    exchange_atoms(i) = C_atoms(1);
end
%
Secretion_fluxes = exchange_fluxes(exchange_fluxes > 0.0);
Secretion_atoms = exchange_atoms(exchange_fluxes > 0.0);
Uptake_fluxes = exchange_fluxes(exchange_fluxes < 0.0);
Uptake_atoms = exchange_atoms(exchange_fluxes < 0.0);
%
Uptake_tmp = abs(sum(Uptake_fluxes.*Uptake_atoms));
Secretion_tmp = abs(sum(Secretion_fluxes.*Secretion_atoms));
%
iEK1008_14_C = (Uptake_tmp - Secretion_tmp)/Uptake_tmp; 



