% load in the file
fileName = 'iYO844.mat';
if ~exist('modelOri','var')
modelOri = readCbModel(fileName);
end

% change all ATP dependent reactions bounds to 1000
modelOri = changeRxnBounds(modelOri,'ATPM',1000,'u');
model = modelOri;

% check if this worked
printConstraints(model)

% set maximum biomass production as the objective function
% need to change this line to identify the biomass reaction in the mat file
model = changeObjective(model,'Biomass_Ecoli_core_N(w/GAM)-Nmet2');


