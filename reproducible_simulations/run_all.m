function run_all(groups)
% RUN_ALL  Re-run the experiments and save every figure they open.
%
%   run_all                         % all groups listed below
%   run_all({'01_uniqueness'})      % only selected groups
%
% Figures go to results/output/<group>/<script>_figN.png (+ .fig in MATLAB),
% console output to results/output/<group>/<script>.log, and the software
% environment to results/output/environment.txt.
% Each script runs inside its own output folder, so files it saves itself
% (saveas / exportgraphics) also land there.

root = setup_paths();

if nargin < 1
    groups = {'01_uniqueness', '02_kingmaker', '03_optimal_control', '04_real_data'};
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
% Separate workspace: a 'clear' inside the script only wipes this function.
cd(outDir);
run(scriptPath);
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
