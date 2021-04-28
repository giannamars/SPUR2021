% clear all removes the solver configuration!
%
initCobraToolbox(false) % false, since don't want to update
% 
changeCobraSolver('gurobi','all');   % change LP solver to GUROBI
%
% load the model
fileName = 'iAF692.mat';
if ~exist('modelOri','var')
    modelOri = readCbModel(fileName);
end
%
model = modelOri; % keep original copy in case we modify the metabolic network
outmodel = writeCbModel(model, 'xls', 'iAF692_model.xls'); % write all reactions to xls file
%
% Set the lower bound of all exchange and sink (siphon) reactions to ensure that only
% those metabolites that are supposed to be taken up are indded supplied to the model.
modelexchanges1 = strmatch('Ex_', model.rxns); 
modelexchanges4 = strmatch('EX_', model.rxns); % identify 70 exchange reactions
modelexchanges2 = strmatch('DM_', model.rxns);
modelexchanges3 = strmatch('sink_', model.rxns);
% To make sure we've found all exchange reactions via string matching, we
% also check for reactions that contain only one non-zero entry in the S matrix column
selExc = (find(full((sum(abs(model.S)==1, 1)==1) & (sum(model.S~=0) == 1))))'; % identify 70 exchange reactions
%
% set lower bounds of exchange reactions to zero
modelexchanges = unique([modelexchanges1; modelexchanges2; modelexchanges3; modelexchanges4; selExc; BM]);
model.lb(find(ismember(model.rxns, model.rxns(modelexchanges))))=0;
%
% set uber bounds to 1000
model.ub(selExc) = 1000; 
%
% biomass optimization
model = changeObjective(model,'BIOMASS_Mb_30');

%%% Test for growth on different sources 
% Glucose under aerobic conditions


model.lb(find(ismember(model.rxns, 'EX_h2o[e]'))) = -1000;
model.ub(find(ismember(model.rxns, 'EX_h2o[e]'))) = 1000;
model.ub(find(ismember(model.rxns, 'EX_co2[e]'))) = 1000;
model.lb(find(ismember(model.rxns, 'EX_glc_D[e]'))) = -1;  % constrain uptake flux
model.ub(find(ismember(model.rxns, 'EX_glc_D[e]'))) = -1;

% solve
FBAsolution = optimizeCbModel(model,'max');

