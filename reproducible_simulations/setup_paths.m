function root = setup_paths()
% SETUP_PATHS  Add src/ to the path; call once before running a single experiment,
%   e.g.  >> setup_paths, cd experiments/fig3_uniqueness, adop_multi_unique

root = fileparts(mfilename('fullpath'));
addpath(fullfile(root, 'src'));
end
