%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           DDPC using an Optimal-IV
%           Authors: R. Dinkla, T. Oomen, J.W. van Wingerden
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
opts = init_opts(N=1e4);

[Re, N, f, Ncl] = deal(opts.Re, opts.N, opts.f, opts.Ncl);

% Requirements:
% 1) Casadi                                     v3.6.7
% 2) Control System Toolbox                     v24.2
% 3) Robust Control Toolbox                     v24.2
% 4) Statistics and Machine Learning Toolbox    v24.2
% 5) crameri_colours                            v1.09
%    Used for plotting if opts.plot = true. Obtained from
%    https://nl.mathworks.com/matlabcentral/fileexchange/68546-crameri-perceptually-uniform-scientific-colormaps
% 6) Parallel Computing Toolbox                 v24.2

% ------------------------- add relevant paths ----------------------------
add_paths(opts);

%% Simulation settings
fprintf('Setting simulation settings...\n');

% ================== initialize simulations ===============================
% get plant, initial controller (Cz0), CL-system (Tcl0), signal names, etc.
[plant,nu,ny,Cz0,Tcl0,opts,sigs] = init_sims(opts);

% ============ set p values to iterate over & number of seeds per p =======
% set p values to iterate over
pmin = max(ss2lag(plant),ss2lag(Cz0)); % take max -> if rho > p approx_IV methods deliver shorter IVs
pmax = 50;
nP   = 2; % number of p values to iterate over
p_all = ceil(linspace(pmin,pmax,nP));

% set seeds to use for iterations
spP = 100;
seeds = reshape(1:nP*spP,spP,nP);

% ================== saving data and settings =============================
% saving this data in data\raw\sys#\dp\<subdir1>
src_dir = pwd;
cd('..'); proj_dir = pwd;
sys_dir = fullfile(pwd,'data','raw',sprintf('sys%d',opts.sys));  % -> data\raw\sys#
if ~isfolder(sys_dir)
    mkdir(sys_dir);
end
cd(src_dir);

% create data\raw\sys#\dp if it doesn't exist yet
if ~isfolder(fullfile(sys_dir,'dp'))
    mkdir(fullfile(sys_dir,'dp'))
end

% create subdir1
subdir1 = name_subdir1(pmin,pmax,nP,opts); % subdir1 name
subdir1 = fullfile(sys_dir,'dp',subdir1);  % subdir1 path
mkdir(subdir1);

% copy dependent .m files to data\raw\sys#\dp\<subdir1>\mfiles
copy_dependencies(src_dir,subdir1,'main_dp.m');

% save overall settings to data\raw\sys#\dp\<subdir1>\dp_settings.mat
save(fullfile(subdir1,'dp_settings.mat'),'pmin','pmax','nP','p_all','spP','seeds','plant','nu','ny','Cz0','Tcl0','opts','sigs');

%% ========================== iterate over p and seeds ====================
if ismember('SlurmProfile1',parallel.clusterProfiles) % on cluster?
    % submit_dp_jobs(subdir1, nP, spP, proj_dir, opts);

    myCluster = parcluster('SlurmProfile1');
    ntasksTotal = nP * spP;
    MaxTasksPerJob = 60;
    nJobs = ceil(ntasksTotal/MaxTasksPerJob);
    idx1 = 1;
    [dTtot, dTmax, ntasks_old] = deal(0);
    nMins = 15;

    % iterate over jobs
    for iJob = 1:nJobs
        tStart = tic;

        % determine number of tasks for job
        if iJob < nJobs
            ntasks = MaxTasksPerJob;
        else
            ntasks = ntasksTotal - (nJobs-1)*MaxTasksPerJob; % however many tasks remain
        end
        idxs = idx1:(idx1+ntasks-1); % ntasks -> task indices

        % estimating time left
        Tleft = nMins - (dTtot + dTmax*1.5)/60; % time left [min] after parfor
        Tleft_str = string(duration(0,Tleft,0,'Format','mm:ss'));
        nMins_str = string(duration(0,nMins,0,'Format','mm:ss'));
        dTtot_str = string(duration(0,0,dTtot,'Format','mm:ss'));
        dTmax_str = string(duration(0,0,dTmax,'Format','mm:ss'));
        
        clc;
        fprintf("Starting job %d\n",iJob);
        fprintf('Number of tasks %d (current), %d (old)\n',ntasks,ntasks_old);
        fprintf('Time left: %s = %s - %s - 1.5*%s\n',Tleft_str,nMins_str,dTtot_str,dTmax_str);

        % reinitialize pool if little time is left or number of tasks changed (possible at last iteration)
        if ntasks ~= ntasks_old || Tleft <= 0
            fprintf('re-initializing pool\n');
            if ~isempty(gcp('nocreate'))
                delete(gcp('nocreate'));
            end
            myCluster.SubmitArguments = SlurmSubmitArgs('dp',nMins,ntasks=ntasks,cpt=1,GB=3.8);
            dTtot  = 0;
            tStart = tic;
            parpool(myCluster,ntasks);
        end

        % iterate over tasks
        parfor iii = 1:ntasks
            idx = idxs(iii)
            run_p(idx,opts,spP,nP,seeds,p_all,f,N,Ncl,ny,nu,Re,plant,subdir1,sigs,Cz0,Tcl0,proj_dir)
        end

        % update timer:
        dT = toc(tStart);       % time taken for this job
        dTtot = dTtot + dT;     % est. of total time used slurm job
        dTmax = max(dTmax,dT);  % est. of max. time used for a parfor iteration
        
        % update others
        idx1 = idxs(end) + 1;
        ntasks_old = ntasks;
    end

else
    fprintf('using the local profile');
    if isempty(gcp('nocreate'))
        myCluster = parcluster('local');
        nworker = myCluster.NumWorkers; % (max.) workers per node
        parpool(myCluster,nworker);
    end
    parfor iii = 1:nP*spP
        run_p(iii,opts,spP,nP,seeds,p_all,f,N,Ncl,ny,nu,Re,plant,subdir1,sigs,Cz0,Tcl0,proj_dir)
    end
end

%% Helper functions
% set name of subdir 1
function subdir1 = name_subdir1(pmin,pmax,nP,opts)
[Re, N] = deal(opts.Re, opts.N);

% Helper function to trim to minimal digits in scientific notation
trimmed_exp = @(x) regexprep(sprintf('%e', x), '(\.\d*?)0+(e[+-]?\d+)', '$1$2'); % trims trailing 0s
% Also remove . if nothing follows
trimmed_exp = @(x) regexprep(trimmed_exp(x), '\.(e)', '$1');

% Apply formatting
N_s   = trimmed_exp(N);
Re_s  = trimmed_exp(Re);
subdir1 = sprintf('p_%d_%d_%d_Re_%s_N_%s',pmin,pmax,nP,Re_s,N_s);
subdir1 = replace(subdir1,'.','p');
subdir1 = replace(subdir1,'+','');
end

function opts = init_opts(opts)
arguments
opts.Re   (1,1) double  = 1e-1;  % innovation noise variance
opts.plot       logical = false;
opts.N    (1,1) double  = 1e3;      % number of Hankel data matrix columns
opts.f    (1,1) double  = 20;
opts.Ncl  (1,1) double  = 1500;     % simulation length of SPC
opts.dRk  (1,1) double  = 1;        % weights
opts.Rk   (1,1) double  = 1;
opts.Qk   (1,1) double  = 1e2;
opts.save       logical = true;     % save data
opts.sys  (1,1) double = 1;         % flag for model selection
end
end

%% Submit Slurm job arrays for DP simulations
function submit_dp_jobs(subdir1, nP, spP, proj_dir, opts)
%SUBMIT_DP_JOBS  Create and submit SLURM job arrays for DP simulations.
%
% Usage:
%   submit_dp_jobs(subdir1, nP, spP, proj_dir, opts)
%
% Inputs:
%   subdir1   - Path to subdirectory for this DP run (string)
%   nP        - Number of p-values (integer)
%   spP       - Seeds per p (integer)
%   proj_dir  - Project root directory (string)
%   opts      - Struct containing fields: f, N, Ncl, Re, sys (for naming)

fprintf('\n=== Submitting %d job arrays to SLURM ===\n', nP);

dp_dir = fullfile(proj_dir,'data','raw',sprintf('sys%d',opts.sys),'dp');  % typically data/raw/sys#/dp

% Ensure dp/ipX directories exist
for iP = 1:nP
    % define and create directories
    ip_dir = fullfile(dp_dir, sprintf('ip%d', iP));
    out_dir = fullfile(ip_dir,'out');
    err_dir = fullfile(ip_dir,'err');
    if ~isfolder(ip_dir)
        mkdir(ip_dir);
    end
    if ~isfolder(out_dir)
        mkdir(out_dir);
    end
    if ~isfolder(err_dir)
        mkdir(err_dir);
    end

    %% --- Create SLURM job array submission script for this iP ---
    script_name = fullfile(ip_dir, sprintf('submit_job_ip%d.sh', iP));
    script_header = sprintf(...
        ['#!/bin/bash\n' ...
        '#SBATCH --job-name=ip%d_s%d\n' ...               iP, opts.sys
        '#SBATCH --partition=compute\n' ...
        '#SBATCH --time=00:15:00\n' ...
        '#SBATCH --ntasks=1\n' ...
        '#SBATCH --cpus-per-task=1\n' ...
        '#SBATCH --mem-per-cpu=3900M\n' ...
        '#SBATCH --account=research-me-dcsc\n' ...
        '#SBATCH --array=1-%d\n' ...                        spP
        '#SBATCH --output=%s/out/job.%%A_%%a.out\n' ...     ip_dir
        '#SBATCH --error=%s/err/job.%%A_%%a.err\n\n'],...   ip_dir
        iP, opts.sys, spP, ip_dir, ip_dir);
    batch_commands = sprintf(...
        ['module load matlab\n\n' ...
        'task_id=$SLURM_ARRAY_TASK_ID\n' ...
        'seed_idx=$(( (%d - 1)*%d + task_id ))\n\n'],...    iP, spP
        iP, spP);
    matlab_commands = sprintf(...
        ['matlab -nosplash -nodesktop -r "'...
        'try; ' ...
            'addpath(genpath(pwd)); '...     should add src & its subdirs to path
            'subdir1=''%s''; '...            define subdir1                                                                                             subdir1
            'addpath(genpath(subdir1)); '... add subdir1 to path
            'fprintf(''loading settings in subdir1\\n''); load(fullfile(subdir1,''dp_settings.mat'')); ' ...        loading settings
            'opts, spP, nP, seeds, p_all, ny, nu, plant, sigs, Cz0, Tcl0, ' ...                                     display used loaded variables
            'fprintf(''defining variables\\n''); ii=$seed_idx, f=%d, N=%d, Ncl=%d, Re=%g, proj_dir=''%s'', ' ...    defining other used variables       opts.f, opts.N, opts.Ncl, opts.Re, proj_dir
            'addpath(genpath(fullfile(proj_dir,''bin''))); '... needed for Casadi
            'run_p(ii,opts,spP,nP,seeds,p_all,f,N,Ncl,ny,nu,Re,plant,subdir1,sigs,Cz0,Tcl0,proj_dir); ' ...         execute run with p
        'catch ME; '...
            'disp(getReport(ME)); exit(1); '...
        'end; '...
        'exit;"\n'], ...
        subdir1, opts.f, opts.N, opts.Ncl, opts.Re, proj_dir);
    script_footer = sprintf('echo "Job ip%d_$task_id completed."', iP);
    script_content = [script_header, batch_commands, matlab_commands, script_footer];

    % Write script file
    if isfile(script_name)
        delete(script_name);
    end
    fid = fopen(script_name, 'w');
    fprintf(fid, '%s', script_content);
    fclose(fid);
    fileattrib(script_name, '+x');

    %% --- Add cleanup script for after all array jobs finish ---
    cleanup_script = fullfile(ip_dir, sprintf('cleanup_ip%d.sh', iP));
    cleanup_content = sprintf([...
    '#!/bin/bash\n',...
    '#SBATCH --job-name=CLip%d_s%d\n',...   iP, opts.sys
    '#SBATCH --partition=compute\n',...
    '#SBATCH --time=00:05:00\n',...
    '#SBATCH --ntasks=1\n',...
    '#SBATCH --cpus-per-task=1\n',...
    '#SBATCH --mem-per-cpu=500M\n',...
    '#SBATCH --output=%s/cleanup.out\n',...         ip_dir
    '#SBATCH --error=%s/cleanup.err\n\n',...        ip_dir
    'ip_dir="%s"\n',...                             ip_dir
    'out_dir="$ip_dir/out"\n',...
    'err_dir="$ip_dir/err"\n',...
    'ip_name=$(basename "${ip_dir}")   # e.g. "ip1"\n',...
    '\n',...
    'echo "Checking job outputs in $out_dir ..."\n',...
    '\n',...
    'for out_file in "$out_dir"/job.*.out; do\n',...
    '  [ -f "$out_file" ] || continue\n',...
    '  filename=$(basename "$out_file")\n',...
    '  job_desc=${filename#job.}\n',...
    '  job_desc=${job_desc%%.out}\n',...
    '  taskid=${job_desc#*_}\n',...
    '\n',...
    '  # Check the last few lines for both success messages\n',...
    '  if tail -n 5 "$out_file" | grep -q "File saved successfully" && \\\n',...
    '     tail -n 5 "$out_file" | grep -q "Job ${ip_name}_${taskid} completed"; then\n',...
    '    rm -f "$out_file"\n',...
    '    rm -f "$err_dir/job.${job_desc}.err"\n',...
    '    echo "Cleaned logs for job ${job_desc}"\n',...
    '  else\n',...
    '    echo "Keeping logs for job ${job_desc} (not confirmed)."\n',...
    '  fi\n',...
    'done\n',...
    '\n',...
    'remaining_out=$(find "$out_dir" -maxdepth 1 -type f | wc -l)\n',...
    'remaining_err=$(find "$err_dir" -maxdepth 1 -type f | wc -l)\n',...
    '\n',...
    'if [ "$remaining_out" -eq 0 ] && [ "$remaining_err" -eq 0 ]; then\n',...
    '  echo "All job logs cleaned. Removing $ip_dir"\n',...
    '  rm -rf "$ip_dir"\n',...
    'else\n',...
    '  echo "Logs remain in $ip_dir; not deleting."\n',...
    'fi'], ...
    iP, opts.sys, ip_dir, ip_dir, ip_dir);

    if isfile(cleanup_script)
        delete(cleanup_script);
    end
    fid2 = fopen(cleanup_script, 'w');
    fprintf(fid2, '%s', cleanup_content);
    fclose(fid2);
    fileattrib(cleanup_script, '+x');
    
    %% Submit jobs
    [status, msg] = system(sprintf('sbatch "%s"', script_name));
    if status == 0
        % fprintf('Submitted job array %d/%d: %s\n', iP, nP, strtrim(msg));
        tokens = regexp(msg, '\d+', 'match');
        jobid = str2double(tokens{1});

        % Submit cleanup dependent on successful completion
        dep_cmd = sprintf('sbatch --dependency=afterok:%d "%s"', jobid, cleanup_script); % <- note dependency
        [~, msg2] = system(dep_cmd);
        tokens2 = regexp(msg2, '\d+', 'match');
        jobid2 = str2double(tokens2{1});

        fprintf('Submitted job array %d/%d (JobID %d) with cleaner upon completion (JobID %d)\n', iP, nP, jobid, jobid2);

    else
        warning('Failed to submit job array %d: %s', iP, msg);
    end
end

fprintf('All %d job arrays submitted.\n', nP);
end
