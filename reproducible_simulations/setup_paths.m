function root = setup_paths()
% SETUP_PATHS  Add the shared simulation engines (src/) to the path.
%   Call once per session before running any experiment script on its own:
%       >> setup_paths
%       >> cd experiments/01_uniqueness, adop_multi_unique

root = fileparts(mfilename('fullpath'));
addpath(fullfile(root, 'src'));
end
