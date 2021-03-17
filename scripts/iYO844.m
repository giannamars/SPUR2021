% load in the file
fileName = 'iYO844.mat';
if ~exist('modelOri','var')
model = readCbModel(fileName);
end

% check if reactions bounds are set to -/+ 1000
% if not, how to change them? doc changeRxnBounds
printConstraints(model)

% set maximum biomass production as the objective function
% need to change this line to identify the biomass reaction in the mat file
%model = changeObjective(model,'Biomass_Ecoli_core_N(w/GAM)-Nmet2');

model.rxnNames



