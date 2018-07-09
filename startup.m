% adding folders to path
addpath(genpath('helper'));
addpath(genpath('functions'));
addpath(genpath('examples'));
addpath(genpath('gaussUnkownMean'));
addpath(genpath('gaussUnknownVar'));
addpath(genpath('classes'));


% creating results directory if not existing
if ~exist('results','dir')
    mkdir('./results/');
end
