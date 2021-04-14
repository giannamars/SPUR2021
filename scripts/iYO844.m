

% load in the file
fileName = 'iYO844.mat';
if ~exist('modelOri','var')
modelOri = readCbModel(fileName);
end

model = modelOri;

model = changeObjective(model,'BIOMASS_BS_10');

FBAsolution = optimizeCbModel(model,'max');

