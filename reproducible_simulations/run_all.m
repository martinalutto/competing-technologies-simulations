function run_all(groups)
% RUN_ALL  Re-run the experiments and save every figure they open.
%   run_all                       % all groups
%   run_all({'fig3_uniqueness'})  % selected groups
% Figures and logs go to results/output/<group>/, environment info to results/output/.

root = setup_paths();

if nargin < 1
    groups = {'fig1_market_share', 'fig3_uniqueness', 'fig4_entry', 'fig5-7_optimal_control', ...
              'extra_kingmaker_analysis', 'extra_problem4_kingmaker_control'};
end

outRoot = fullfile(root, 'results', 'output');
if ~exist(outRoot, 'dir'), mkdir(outRoot); end
write_environment(fullfile(outRoot, 'environment.txt'));

summary = {};
for g = 1:numel(groups)
    expDir = fullfile(root, 'experiments', groups{g});
    outDir = fullfile(outRoot, groups{g});
    if ~exist(outDir, 'dir'), mkdir(outDir); end

    files = dir(fullfile(expDir, '*.m'));
    for f = 1:numel(files)
        [~, name] = fileparts(files(f).name);
        fprintf('\n>>> %s/%s\n', groups{g}, name);
        [ok, msg, secs] = run_one(fullfile(expDir, files(f).name), outDir, name);
        summary(end+1, :) = {groups{g}, name, ok, secs, msg}; %#ok<AGROW>
    end
end

fprintf('\n================ SUMMARY ================\n');
for r = 1:size(summary, 1)
    status = 'OK';
    if ~summary{r, 3}, status = ['FAILED: ' summary{r, 5}]; end
    fprintf('%-20s %-50s %7.1f s  %s\n', summary{r, 1}, summary{r, 2}, summary{r, 4}, status);
end
end

function [ok, msg, secs] = run_one(scriptPath, outDir, name)
here = pwd;
close all;
logFile = fullfile(outDir, [name '.log']);
if exist(logFile, 'file'), delete(logFile); end
diary(logFile);
ok = true; msg = '';
t0 = tic;
try
    exec_script(scriptPath, outDir);
catch err
    ok = false; msg = err.message;
    fprintf(2, 'ERROR in %s: %s\n', name, msg);
end
secs = toc(t0);
diary off;
cd(here);

figs = findall(0, 'Type', 'figure');
if ~isempty(figs)
    [~, order] = sort(arrayfun(@(h) double(get(h, 'Number')), figs));
    figs = figs(order);
end
for k = 1:numel(figs)
    base = fullfile(outDir, sprintf('%s_fig%d', name, k));
    print(figs(k), [base '.png'], '-dpng', '-r200');
    if ~is_octave(), savefig(figs(k), [base '.fig']); end
end
close all;
end

function exec_script(scriptPath, outDir)
% Own workspace, so a 'clear' in the script is harmless; run() is avoided
% because it changes into the script's folder.
[scriptDir, scriptName] = fileparts(scriptPath);
addpath(scriptDir);
cleanup = onCleanup(@() rmpath(scriptDir));
cd(outDir);
eval(scriptName);
end

function write_environment(file)
fid = fopen(file, 'w');
fprintf(fid, 'Date: %s\n', datestr(now));
if is_octave()
    fprintf(fid, 'GNU Octave %s\n', OCTAVE_VERSION);
else
    fprintf(fid, 'MATLAB %s\n', version);
    v = ver;
    for k = 1:numel(v), fprintf(fid, '  %s %s\n', v(k).Name, v(k).Version); end
end
fprintf(fid, 'Computer: %s\n', computer);
fclose(fid);
end

function tf = is_octave()
tf = exist('OCTAVE_VERSION', 'builtin') ~= 0;
end
