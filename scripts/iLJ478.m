clear;
% load the model
fileName = 'iLJ478.mat';
if ~exist('modelOri','var')
    modelOri = readCbModel(fileName);
end
model = modelOri; % keep original copy in case we modify the metabolic network
outmodel = writeCbModel(model, 'xls', 'iLJ478_model.xls'); % write all reactions to xls file
modelClosed = model;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%modelClosed = model;
% Find all exchange reactions in the model
modelexchanges1 = strmatch('Ex_', modelClosed.rxns); 
modelexchanges4 = strmatch('EX_', modelClosed.rxns);          % identify 70 exchange reactions
modelexchanges2 = strmatch('DM_', modelClosed.rxns);
modelexchanges3 = strmatch('sink_', modelClosed.rxns);
% To make sure we've found all exchange reactions via string matching, we
% also check for reactions that contain only one non-zero entry in the S matrix column
selExc = (find(full((sum(abs(modelClosed.S)==1, 1)==1) & (sum(modelClosed.S~=0) == 1))))'; % identify 70 exchange reactions
modelexchanges = unique([modelexchanges1; modelexchanges2; modelexchanges3; modelexchanges4; selExc]);
reaction_abbrvs = modelClosed.rxns(modelexchanges); % store reaction abbreviations in model
reaction_names = modelClosed.rxnNames(modelexchanges); % store reaction names
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For the set of 13 models, we want to manipulate the availability of 14
% individual C sources in an aerobic environment, and calculate CUE under exclusive uptake of each
% metabolite separately
%
% set lower bounds (=uptake) of all exchange reactions to zero, meaning to uptake
modelClosed.lb(find(ismember(modelClosed.rxns, modelClosed.rxns(modelexchanges))))=0;
% set uber bounds(=excretion) to 1000 just to be safe
%modelClosed.ub(selExc) = 1000; 
%modelClosed.lb(selExc) = -100;
% change objective function to biomass optimization
modelClosed = changeObjective(modelClosed,'BIOMASS_Ecoli_TM');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Acetate
modelClosedOri = modelClosed;
modelClosed.lb(find(ismember(modelClosed.rxns, 'EX_o2_e'))) = -100;
modelClosed.lb(find(ismember(modelClosed.rxns, 'EX_h2o_e'))) = -1000;
modelClosed.ub(find(ismember(modelClosed.rxns, 'EX_h2o_e'))) = 1000;
modelClosed.lb(find(ismember(modelClosed.rxns, 'EX_co2_e'))) = -1000;
modelClosed.ub(find(ismember(modelClosed.rxns, 'EX_co2_e'))) = 1000;
modelClosed.lb(find(ismember(modelClosed.rxns, 'EX_ac_e'))) = -10;
FBA = optimizeCbModel(modelClosed, 'max');

if isnan(FBA.f)
    iLJ478_1_C = NaN;
else
    exchange_fluxes = FBA.x(find(ismember(model.rxns, reaction_abbrvs)));
    %
    exchange_atoms = zeros(length(exchange_fluxes),1);
    for i = 1:length(exchange_atoms)
        tmp1 = printRxnFormula(modelClosed,reaction_abbrvs(i));
        tmp2 = strtrim(regexp(tmp1{1}, '->|<=>', 'split'));
        metID = findMetIDs(modelClosed, tmp2{1});
        metFormula= modelClosed.metFormulas(metID);
        try
            atomicStructure = parse_formula(metFormula{1});
            C_atoms = zeros(1);
            if(isfield(atomicStructure, 'C'))
                C_atoms(1) = atomicStructure.C;
            else
                C_atoms(1) = 0;
            end
            exchange_atoms(i) = C_atoms(1);
        catch
            warning('Problem using function.  Assigning a value of 0.');
            exchange_atoms(i) = 0.0;
        end
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
    iLJ478_1_C = (Uptake_tmp - Secretion_tmp)/Uptake_tmp;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Pyruvate
modelClosedOri = modelClosed;
modelClosed.lb(find(ismember(modelClosed.rxns, 'EX_o2_e'))) = -1000;
modelClosed.lb(find(ismember(modelClosed.rxns, 'EX_h2o_e'))) = -1000;
modelClosed.ub(find(ismember(modelClosed.rxns, 'EX_h2o_e'))) = 1000;
modelClosed.lb(find(ismember(modelClosed.rxns, 'EX_co2_e'))) = -1000;
modelClosed.ub(find(ismember(modelClosed.rxns, 'EX_co2_e'))) = 1000;
modelClosed.lb(find(ismember(modelClosed.rxns, 'EX_pyr_e'))) = -1;
FBA = optimizeCbModel(modelClosed, 'max');

if isnan(FBA.f)
    iLJ478_2_C = NaN;
else
    exchange_fluxes = FBA.x(find(ismember(modelClosed.rxns, reaction_abbrvs)));%
    %
    exchange_atoms = zeros(length(exchange_fluxes),1);
    for i = 1:length(exchange_atoms)
        tmp1 = printRxnFormula(modelClosed,reaction_abbrvs(i));
        tmp2 = strtrim(regexp(tmp1{1}, '->|<=>', 'split'));
        metID = findMetIDs(modelClosed, tmp2{1});
        metFormula= modelClosed.metFormulas(metID);
        try
            atomicStructure = parse_formula(metFormula{1});
            C_atoms = zeros(1);
            if(isfield(atomicStructure, 'C'))
                C_atoms(1) = atomicStructure.C;
            else
                C_atoms(1) = 0;
            end
            exchange_atoms(i) = C_atoms(1);
        catch
            warning('Problem using function.  Assigning a value of 0.');
            exchange_atoms(i) = 0.0;
        end
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
    iLJ478_2_C = (Uptake_tmp - Secretion_tmp)/Uptake_tmp;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% D-glucose
modelClosedOri = modelClosed;
modelClosed.lb(find(ismember(modelClosed.rxns, 'EX_o2_e'))) = -1000;
modelClosed.lb(find(ismember(modelClosed.rxns, 'EX_h2o_e'))) = -1000;
modelClosed.ub(find(ismember(modelClosed.rxns, 'EX_h2o_e'))) = 1000;
modelClosed.lb(find(ismember(modelClosed.rxns, 'EX_co2_e'))) = -1000;
modelClosed.ub(find(ismember(modelClosed.rxns, 'EX_co2_e'))) = 1000;
modelClosed.lb(find(ismember(modelClosed.rxns, 'EX_glc__D_e'))) = -1;
FBA = optimizeCbModel(modelClosed, 'max');
if isnan(FBA.f)
    iLJ478_3_C = NaN;
else
    exchange_fluxes = FBA.x(find(ismember(modelClosed.rxns, reaction_abbrvs)));%
    %
    exchange_atoms = zeros(length(exchange_fluxes),1);
    for i = 1:length(exchange_atoms)
        tmp1 = printRxnFormula(modelClosed,reaction_abbrvs(i));
        tmp2 = strtrim(regexp(tmp1{1}, '->|<=>', 'split'));
        metID = findMetIDs(modelClosed, tmp2{1});
        metFormula= modelClosed.metFormulas(metID);
            try
                atomicStructure = parse_formula(metFormula{1});
                C_atoms = zeros(1);
                if(isfield(atomicStructure, 'C'))
                    C_atoms(1) = atomicStructure.C;
                else
                    C_atoms(1) = 0;
                end
                exchange_atoms(i) = C_atoms(1);
            catch
                warning('Problem using function.  Assigning a value of 0.');
                exchange_atoms(i) = 0.0;
            end
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
    iLJ478_3_C = (Uptake_tmp - Secretion_tmp)/Uptake_tmp; 
end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fumarate
modelClosedOri = modelClosed;
modelClosed.lb(find(ismember(modelClosed.rxns, 'EX_o2_e'))) = -1000;
modelClosed.lb(find(ismember(modelClosed.rxns, 'EX_h2o_e'))) = -1000;
modelClosed.ub(find(ismember(modelClosed.rxns, 'EX_h2o_e'))) = 1000;
modelClosed.lb(find(ismember(modelClosed.rxns, 'EX_co2_e'))) = -1000;
modelClosed.ub(find(ismember(modelClosed.rxns, 'EX_co2_e'))) = 1000;
modelClosed.lb(find(ismember(modelClosed.rxns, 'EX_fum_e'))) = -1;
FBA = optimizeCbModel(modelClosed, 'max');

if isnan(FBA.f)
    iLJ478_4_C = NaN;
else

    exchange_fluxes = FBA.x(find(ismember(modelClosed.rxns, reaction_abbrvs)));%
    %
    exchange_atoms = zeros(length(exchange_fluxes),1);
    for i = 1:length(exchange_atoms)
        tmp1 = printRxnFormula(modelClosed,reaction_abbrvs(i));
        tmp2 = strtrim(regexp(tmp1{1}, '->|<=>', 'split'));
        metID = findMetIDs(modelClosed, tmp2{1});
        metFormula= modelClosed.metFormulas(metID);

            try
                atomicStructure = parse_formula(metFormula{1});
                C_atoms = zeros(1);
                if(isfield(atomicStructure, 'C'))
                    C_atoms(1) = atomicStructure.C;
                else
                    C_atoms(1) = 0;
                end
                exchange_atoms(i) = C_atoms(1);
            catch
                warning('Problem using function.  Assigning a value of 0.');
                exchange_atoms(i) = 0.0;
            end
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
    iLJ478_4_C = (Uptake_tmp - Secretion_tmp)/Uptake_tmp;  
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Acetaldehyde
modelClosedOri = modelClosed;
modelClosed.lb(find(ismember(modelClosed.rxns, 'EX_o2_e'))) = -1000;
modelClosed.lb(find(ismember(modelClosed.rxns, 'EX_h2o_e'))) = -1000;
modelClosed.ub(find(ismember(modelClosed.rxns, 'EX_h2o_e'))) = 1000;
modelClosed.lb(find(ismember(modelClosed.rxns, 'EX_co2_e'))) = -1000;
modelClosed.ub(find(ismember(modelClosed.rxns, 'EX_co2_e'))) = 1000;
modelClosed.lb(find(ismember(modelClosed.rxns, 'EX_acald_e'))) = -1;
modelClosed.ub(find(ismember(modelClosed.rxns, 'EX_acald_e'))) = -1;
FBA = optimizeCbModel(modelClosed, 'max');

if isnan(FBA.f)
    iLJ478_5_C = NaN;
else
    exchange_fluxes = FBA.x(find(ismember(modelClosed.rxns, reaction_abbrvs)));%
    %
    exchange_atoms = zeros(length(exchange_fluxes),1);
    for i = 1:length(exchange_atoms)
        tmp1 = printRxnFormula(modelClosed,reaction_abbrvs(i));
        tmp2 = strtrim(regexp(tmp1{1}, '->|<=>', 'split'));
        metID = findMetIDs(modelClosed, tmp2{1});
        metFormula= modelClosed.metFormulas(metID);
            try
                atomicStructure = parse_formula(metFormula{1});
                C_atoms = zeros(1);
                if(isfield(atomicStructure, 'C'))
                    C_atoms(1) = atomicStructure.C;
                else
                    C_atoms(1) = 0;
                end
                exchange_atoms(i) = C_atoms(1);
            catch
                warning('Problem using function.  Assigning a value of 0.');
                exchange_atoms(i) = 0.0;
            end
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
    iLJ478_5_C = (Uptake_tmp - Secretion_tmp)/Uptake_tmp;   
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2-oxoglutarate
modelClosedOri = modelClosed;
modelClosed.lb(find(ismember(modelClosed.rxns, 'EX_o2_e'))) = -1000;
modelClosed.lb(find(ismember(modelClosed.rxns, 'EX_h2o_e'))) = -1000;
modelClosed.ub(find(ismember(modelClosed.rxns, 'EX_h2o_e'))) = 1000;
modelClosed.lb(find(ismember(modelClosed.rxns, 'EX_co2_e'))) = -1000;
modelClosed.ub(find(ismember(modelClosed.rxns, 'EX_co2_e'))) = 1000;
modelClosed.lb(find(ismember(modelClosed.rxns, 'EX_akg_e'))) = -1;
modelClosed.ub(find(ismember(modelClosed.rxns, 'EX_akg_e'))) = -1;
FBA = optimizeCbModel(modelClosed, 'max');
if isnan(FBA.f)
    iLJ478_6_C = NaN;
else
exchange_fluxes = FBA.x(find(ismember(modelClosed.rxns, reaction_abbrvs)));%
%
exchange_atoms = zeros(length(exchange_fluxes),1);
for i = 1:length(exchange_atoms)
    tmp1 = printRxnFormula(modelClosed,reaction_abbrvs(i));
    tmp2 = strtrim(regexp(tmp1{1}, '->|<=>', 'split'));
    metID = findMetIDs(modelClosed, tmp2{1});
    metFormula= modelClosed.metFormulas(metID);
            try
                atomicStructure = parse_formula(metFormula{1});
                C_atoms = zeros(1);
                if(isfield(atomicStructure, 'C'))
                    C_atoms(1) = atomicStructure.C;
                else
                    C_atoms(1) = 0;
                end
                exchange_atoms(i) = C_atoms(1);
            catch
                warning('Problem using function.  Assigning a value of 0.');
                exchange_atoms(i) = 0.0;
            end
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
iLJ478_6_C = (Uptake_tmp - Secretion_tmp)/Uptake_tmp;   
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ethanol
modelClosedOri = modelClosed;
modelClosed.lb(find(ismember(modelClosed.rxns, 'EX_o2_e'))) = -1000;
modelClosed.lb(find(ismember(modelClosed.rxns, 'EX_h2o_e'))) = -1000;
modelClosed.ub(find(ismember(modelClosed.rxns, 'EX_h2o_e'))) = 1000;
modelClosed.lb(find(ismember(modelClosed.rxns, 'EX_co2_e'))) = -1000;
modelClosed.ub(find(ismember(modelClosed.rxns, 'EX_co2_e'))) = 1000;
modelClosed.lb(find(ismember(modelClosed.rxns, 'EX_etoh_e'))) = -1;
modelClosed.ub(find(ismember(modelClosed.rxns, 'EX_etoh_e'))) = -1;
FBA = optimizeCbModel(modelClosed, 'max');
if isnan(FBA.f)
    iLJ478_7_C = NaN;
else
exchange_fluxes = FBA.x(find(ismember(modelClosed.rxns, reaction_abbrvs)));%
%
exchange_atoms = zeros(length(exchange_fluxes),1);
for i = 1:length(exchange_atoms)
    tmp1 = printRxnFormula(modelClosed,reaction_abbrvs(i));
    tmp2 = strtrim(regexp(tmp1{1}, '->|<=>', 'split'));
    metID = findMetIDs(modelClosed, tmp2{1});
    metFormula= modelClosed.metFormulas(metID);
            try
                atomicStructure = parse_formula(metFormula{1});
                C_atoms = zeros(1);
                if(isfield(atomicStructure, 'C'))
                    C_atoms(1) = atomicStructure.C;
                else
                    C_atoms(1) = 0;
                end
                exchange_atoms(i) = C_atoms(1);
            catch
                warning('Problem using function.  Assigning a value of 0.');
                exchange_atoms(i) = 0.0;
            end
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
iLJ478_7_C = (Uptake_tmp - Secretion_tmp)/Uptake_tmp;  
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Formate
modelClosedOri = modelClosed;
modelClosed.lb(find(ismember(modelClosed.rxns, 'EX_o2_e'))) = -1000;
modelClosed.lb(find(ismember(modelClosed.rxns, 'EX_h2o_e'))) = -1000;
modelClosed.ub(find(ismember(modelClosed.rxns, 'EX_h2o_e'))) = 1000;
modelClosed.lb(find(ismember(modelClosed.rxns, 'EX_co2_e'))) = -1000;
modelClosed.ub(find(ismember(modelClosed.rxns, 'EX_co2_e'))) = 1000;
modelClosed.lb(find(ismember(modelClosed.rxns, 'EX_for_e'))) = -1;
modelClosed.ub(find(ismember(modelClosed.rxns, 'EX_for_e'))) = -1;
FBA = optimizeCbModel(modelClosed, 'max');
if isnan(FBA.f)
    iLJ478_8_C = NaN;
else
exchange_fluxes = FBA.x(find(ismember(modelClosed.rxns, reaction_abbrvs)));%
%
exchange_atoms = zeros(length(exchange_fluxes),1);
for i = 1:length(exchange_atoms)
    tmp1 = printRxnFormula(modelClosed,reaction_abbrvs(i));
    tmp2 = strtrim(regexp(tmp1{1}, '->|<=>', 'split'));
    metID = findMetIDs(modelClosed, tmp2{1});
            try
                atomicStructure = parse_formula(metFormula{1});
                C_atoms = zeros(1);
                if(isfield(atomicStructure, 'C'))
                    C_atoms(1) = atomicStructure.C;
                else
                    C_atoms(1) = 0;
                end
                exchange_atoms(i) = C_atoms(1);
            catch
                warning('Problem using function.  Assigning a value of 0.');
                exchange_atoms(i) = 0.0;
            end
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
iLJ478_8_C = (Uptake_tmp - Secretion_tmp)/Uptake_tmp;  
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% D-fructose
modelClosedOri = modelClosed;
modelClosed.lb(find(ismember(modelClosed.rxns, 'EX_o2_e'))) = -1000;
modelClosed.lb(find(ismember(modelClosed.rxns, 'EX_h2o_e'))) = -1000;
modelClosed.ub(find(ismember(modelClosed.rxns, 'EX_h2o_e'))) = 1000;
modelClosed.lb(find(ismember(modelClosed.rxns, 'EX_co2_e'))) = -1000;
modelClosed.ub(find(ismember(modelClosed.rxns, 'EX_co2_e'))) = 1000;
modelClosed.lb(find(ismember(modelClosed.rxns, 'EX_fru_e'))) = -1;
modelClosed.ub(find(ismember(modelClosed.rxns, 'EX_fru_e'))) = -1;
FBA = optimizeCbModel(modelClosed, 'max');
if isnan(FBA.f)
    iLJ478_9_C = NaN;
else
exchange_fluxes = FBA.x(find(ismember(modelClosed.rxns, reaction_abbrvs)));%
%
exchange_atoms = zeros(length(exchange_fluxes),1);
for i = 1:length(exchange_atoms)
    tmp1 = printRxnFormula(modelClosed,reaction_abbrvs(i));
    tmp2 = strtrim(regexp(tmp1{1}, '->|<=>', 'split'));
    metID = findMetIDs(modelClosed, tmp2{1});
    metFormula= modelClosed.metFormulas(metID);
            try
                atomicStructure = parse_formula(metFormula{1});
                C_atoms = zeros(1);
                if(isfield(atomicStructure, 'C'))
                    C_atoms(1) = atomicStructure.C;
                else
                    C_atoms(1) = 0;
                end
                exchange_atoms(i) = C_atoms(1);
            catch
                warning('Problem using function.  Assigning a value of 0.');
                exchange_atoms(i) = 0.0;
            end
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
iLJ478_8_C = (Uptake_tmp - Secretion_tmp)/Uptake_tmp;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% L-glutamine
modelClosedOri = modelClosed;
modelClosed.lb(find(ismember(modelClosed.rxns, 'EX_o2_e'))) = -1000;
modelClosed.lb(find(ismember(modelClosed.rxns, 'EX_h2o_e'))) = -1000;
modelClosed.ub(find(ismember(modelClosed.rxns, 'EX_h2o_e'))) = 1000;
modelClosed.lb(find(ismember(modelClosed.rxns, 'EX_co2_e'))) = -1000;
modelClosed.ub(find(ismember(modelClosed.rxns, 'EX_co2_e'))) = 1000;
modelClosed.lb(find(ismember(modelClosed.rxns, 'EX_gln__L_e'))) = -1;
modelClosed.ub(find(ismember(modelClosed.rxns, 'EX_gln__L_e'))) = -1;
FBA = optimizeCbModel(modelClosed, 'max');
if isnan(FBA.f)
    iLJ478_19_C = NaN;
else
exchange_fluxes = FBA.x(find(ismember(modelClosed.rxns, reaction_abbrvs)));%
%
exchange_atoms = zeros(length(exchange_fluxes),1);
for i = 1:length(exchange_atoms)
    tmp1 = printRxnFormula(modelClosed,reaction_abbrvs(i));
    tmp2 = strtrim(regexp(tmp1{1}, '->|<=>', 'split'));
    metID = findMetIDs(modelClosed, tmp2{1});
    metFormula= modelClosed.metFormulas(metID);
            try
                atomicStructure = parse_formula(metFormula{1});
                C_atoms = zeros(1);
                if(isfield(atomicStructure, 'C'))
                    C_atoms(1) = atomicStructure.C;
                else
                    C_atoms(1) = 0;
                end
                exchange_atoms(i) = C_atoms(1);
            catch
                warning('Problem using function.  Assigning a value of 0.');
                exchange_atoms(i) = 0.0;
            end
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
iLJ478_10_C = (Uptake_tmp - Secretion_tmp)/Uptake_tmp;  
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% L-glutamate
modelClosedOri = modelClosed;
modelClosed.lb(find(ismember(modelClosed.rxns, 'EX_o2_e'))) = -1000;
modelClosed.lb(find(ismember(modelClosed.rxns, 'EX_h2o_e'))) = -1000;
modelClosed.ub(find(ismember(modelClosed.rxns, 'EX_h2o_e'))) = 1000;
modelClosed.lb(find(ismember(modelClosed.rxns, 'EX_co2_e'))) = -1000;
modelClosed.ub(find(ismember(modelClosed.rxns, 'EX_co2_e'))) = 1000;
modelClosed.lb(find(ismember(modelClosed.rxns, 'EX_glu__L_e'))) = -1;
modelClosed.ub(find(ismember(modelClosed.rxns, 'EX_glu__L_e'))) = -1;
FBA = optimizeCbModel(modelClosed, 'max');
if isnan(FBA.f)
    iLJ478_11_C = NaN;
else
exchange_fluxes = FBA.x(find(ismember(modelClosed.rxns, reaction_abbrvs)));%
%
exchange_atoms = zeros(length(exchange_fluxes),1);
for i = 1:length(exchange_atoms)
    tmp1 = printRxnFormula(modelClosed,reaction_abbrvs(i));
    tmp2 = strtrim(regexp(tmp1{1}, '->|<=>', 'split'));
    metID = findMetIDs(modelClosed, tmp2{1});
    metFormula= modelClosed.metFormulas(metID);
            try
                atomicStructure = parse_formula(metFormula{1});
                C_atoms = zeros(1);
                if(isfield(atomicStructure, 'C'))
                    C_atoms(1) = atomicStructure.C;
                else
                    C_atoms(1) = 0;
                end
                exchange_atoms(i) = C_atoms(1);
            catch
                warning('Problem using function.  Assigning a value of 0.');
                exchange_atoms(i) = 0.0;
            end
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
iLJ478_11_C = (Uptake_tmp - Secretion_tmp)/Uptake_tmp; 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% D-lactate
modelClosedOri = modelClosed;
modelClosed.lb(find(ismember(modelClosed.rxns, 'EX_o2_e'))) = -1000;
modelClosed.lb(find(ismember(modelClosed.rxns, 'EX_h2o_e'))) = -1000;
modelClosed.ub(find(ismember(modelClosed.rxns, 'EX_h2o_e'))) = 1000;
modelClosed.lb(find(ismember(modelClosed.rxns, 'EX_co2_e'))) = -1000;
modelClosed.ub(find(ismember(modelClosed.rxns, 'EX_co2_e'))) = 1000;
modelClosed.lb(find(ismember(modelClosed.rxns, 'EX_lac__D_e'))) = -1;
modelClosed.ub(find(ismember(modelClosed.rxns, 'EX_lac__D_e'))) = -1;
FBA = optimizeCbModel(modelClosed, 'max');
if isnan(FBA.f)
    iLJ478_12_C = NaN;
else
exchange_fluxes = FBA.x(find(ismember(modelClosed.rxns, reaction_abbrvs)));%
%
exchange_atoms = zeros(length(exchange_fluxes),1);
for i = 1:length(exchange_atoms)
    tmp1 = printRxnFormula(modelClosed,reaction_abbrvs(i));
    tmp2 = strtrim(regexp(tmp1{1}, '->|<=>', 'split'));
    metID = findMetIDs(modelClosed, tmp2{1});
    metFormula= modelClosed.metFormulas(metID);
            try
                atomicStructure = parse_formula(metFormula{1});
                C_atoms = zeros(1);
                if(isfield(atomicStructure, 'C'))
                    C_atoms(1) = atomicStructure.C;
                else
                    C_atoms(1) = 0;
                end
                exchange_atoms(i) = C_atoms(1);
            catch
                warning('Problem using function.  Assigning a value of 0.');
                exchange_atoms(i) = 0.0;
            end
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
iLJ478_12_C = (Uptake_tmp - Secretion_tmp)/Uptake_tmp; 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% L-malate
modelClosedOri = modelClosed;
modelClosed.lb(find(ismember(modelClosed.rxns, 'EX_o2_e'))) = -1000;
modelClosed.lb(find(ismember(modelClosed.rxns, 'EX_h2o_e'))) = -1000;
modelClosed.ub(find(ismember(modelClosed.rxns, 'EX_h2o_e'))) = 1000;
modelClosed.lb(find(ismember(modelClosed.rxns, 'EX_co2_e'))) = -1000;
modelClosed.ub(find(ismember(modelClosed.rxns, 'EX_co2_e'))) = 1000;
modelClosed.lb(find(ismember(modelClosed.rxns, 'EX_mal__L_e'))) = -1;
modelClosed.ub(find(ismember(modelClosed.rxns, 'EX_mal__L_e'))) = -1;
FBA = optimizeCbModel(modelClosed, 'max');
if isnan(FBA.f)
    iLJ478_13_C = NaN;
else
exchange_fluxes = FBA.x(find(ismember(modelClosed.rxns, reaction_abbrvs)));%
%
exchange_atoms = zeros(length(exchange_fluxes),1);
for i = 1:length(exchange_atoms)
    tmp1 = printRxnFormula(modelClosed,reaction_abbrvs(i));
    tmp2 = strtrim(regexp(tmp1{1}, '->|<=>', 'split'));
    metID = findMetIDs(modelClosed, tmp2{1});
    metFormula= modelClosed.metFormulas(metID);
            try
                atomicStructure = parse_formula(metFormula{1});
                C_atoms = zeros(1);
                if(isfield(atomicStructure, 'C'))
                    C_atoms(1) = atomicStructure.C;
                else
                    C_atoms(1) = 0;
                end
                exchange_atoms(i) = C_atoms(1);
            catch
                warning('Problem using function.  Assigning a value of 0.');
                exchange_atoms(i) = 0.0;
            end
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
iLJ478_13_C = (Uptake_tmp - Secretion_tmp)/Uptake_tmp; 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Succinate
modelClosedOri = modelClosed;
modelClosed.lb(find(ismember(modelClosed.rxns, 'EX_o2_e'))) = -1000;
modelClosed.lb(find(ismember(modelClosed.rxns, 'EX_h2o_e'))) = -1000;
modelClosed.ub(find(ismember(modelClosed.rxns, 'EX_h2o_e'))) = 1000;
modelClosed.lb(find(ismember(modelClosed.rxns, 'EX_co2_e'))) = -1000;
modelClosed.ub(find(ismember(modelClosed.rxns, 'EX_co2_e'))) = 1000;
modelClosed.lb(find(ismember(modelClosed.rxns, 'EX_succ_e'))) = -1;
modelClosed.ub(find(ismember(modelClosed.rxns, 'EX_succ_e'))) = -1;
FBA = optimizeCbModel(modelClosed, 'max');
if isnan(FBA.f)
    iLJ478_14_C = NaN;
else
exchange_fluxes = FBA.x(find(ismember(modelClosed.rxns, reaction_abbrvs)));%
%
exchange_atoms = zeros(length(exchange_fluxes),1);
for i = 1:length(exchange_atoms)
    tmp1 = printRxnFormula(modelClosed,reaction_abbrvs(i));
    tmp2 = strtrim(regexp(tmp1{1}, '->|<=>', 'split'));
    metID = findMetIDs(modelClosed, tmp2{1});
    metFormula= modelClosed.metFormulas(metID);
            try
                atomicStructure = parse_formula(metFormula{1});
                C_atoms = zeros(1);
                if(isfield(atomicStructure, 'C'))
                    C_atoms(1) = atomicStructure.C;
                else
                    C_atoms(1) = 0;
                end
                exchange_atoms(i) = C_atoms(1);
            catch
                warning('Problem using function.  Assigning a value of 0.');
                exchange_atoms(i) = 0.0;
            end
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
iLJ478_14_C = (Uptake_tmp - Secretion_tmp)/Uptake_tmp; 
end



