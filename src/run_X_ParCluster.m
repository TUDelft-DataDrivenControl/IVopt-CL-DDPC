function run_X_ParCluster(opts,vs,data_type, optvar)
arguments
    opts struct
    vs struct
    data_type char {mustBeMember(data_type,{'p','Re','N'})} % which variable to iterate over
    optvar.MaxTasksPerJob (1,1) double {mustBeFinite,mustBeReal,mustBePositive, mustBeInteger} = 50 % max. tasks per slurm job
    optvar.nMins (1,1) double {mustBeFinite,mustBeReal,mustBePositive} = 15 % time limit per slurm job [min]
    optvar.GB (1,1) double {mustBeFinite,mustBeReal,mustBePositive} = 3.8 % memory per worker [GB]
end

% unpacking input variables depending on data_type
switch data_type
    case 'p'
               [spP,    nP,    seeds,    p_all,    f,    N,    Ncl,    ny,    nu,    Re,    plant,    subdir1,    sigs,    Cz0,    Tcl0,    proj_dir] = ...
        deal(vs.spP, vs.nP, vs.seeds, vs.p_all, vs.f, vs.N, vs.Ncl, vs.ny, vs.nu, vs.Re, vs.plant, vs.subdir1, vs.sigs, vs.Cz0, vs.Tcl0, vs.proj_dir);
        nX = nP; spX = spP;
    case 'N'
               [spN,    nN,    seeds,    N_all,    p,    f,    Ncl,    ny,    nu,    Re,    plant,    subdir1,    sigs,    Cz0,    Tcl0,    proj_dir] = ...
        deal(vs.spN, vs.nN, vs.seeds, vs.N_all, vs.p, vs.f, vs.Ncl, vs.ny, vs.nu, vs.Re, vs.plant, vs.subdir1, vs.sigs, vs.Cz0, vs.Tcl0, vs.proj_dir);
        nX = nN; spX = spN;
    case 'Re'
               [spRe,    nRe,    seeds,    Re_all,    Ncl,    Nbar,    ny,    plant,    subdir1,    sigs,    Cz0,    Tcl0,    yr0,    yr1,    ur1,    proj_dir] = ...
        deal(vs.spRe, vs.nRe, vs.seeds, vs.Re_all, vs.Ncl, vs.Nbar, vs.ny, vs.plant, vs.subdir1, vs.sigs, vs.Cz0, vs.Tcl0, vs.yr0, vs.yr1, vs.ur1, vs.proj_dir);
        nX = nRe; spX = spRe;
end
[MaxTasksPerJob, nMins, GB] = deal(optvar.MaxTasksPerJob, optvar.nMins, optvar.GB);

%% iterating over jobs
myCluster = parcluster('SlurmProfile1');
ntasksTotal = nX * spX;
nJobs = ceil(ntasksTotal/MaxTasksPerJob);
idx1 = 1;
[dTtot, dTmax, ntasks_old,iSlurmJob] = deal(0);

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
        fprintf('(re-)initializing pool\n');
        if ~isempty(gcp('nocreate'))
            delete(gcp('nocreate'));
        end
        iSlurmJob = iSlurmJob + 1; % counter for slurm job submission
        jobName = sprintf('d%s_%d',data_type,iSlurmJob);
        myCluster.SubmitArguments = SlurmSubmitArgs(jobName,nMins,ntasks=ntasks,cpt=1,GB=GB);
        dTtot  = 0;
        tStart = tic; % restart timer
        parpool(myCluster,ntasks);
    end

    % iterate over tasks
    switch data_type
        case 'p'
            parfor iii = 1:ntasks
                idx = idxs(iii);
                run_p(idx,opts,spP,nP,seeds,p_all,f,N,Ncl,ny,nu,Re,plant,subdir1,sigs,Cz0,Tcl0,proj_dir);
            end
        case 'N'
            parfor iii = 1:ntasks
                idx = idxs(iii);
                run_N(idx,opts,spN,nN,seeds,N_all,p,f,Ncl,ny,nu,Re,plant,subdir1,sigs,Cz0,Tcl0,proj_dir);
            end
        case 'Re'
            parfor iii = 1:ntasks
                idx = idxs(iii);
                run_Re(idx,opts,spRe,nRe,seeds,Re_all,Ncl,Nbar,ny,plant,subdir1,sigs,Cz0,Tcl0,yr0,yr1,ur1,proj_dir);
            end
    end

    % update timer:
    dT = toc(tStart);       % time taken for this job
    dTtot = dTtot + dT;     % est. of total time used slurm job
    dTmax = max(dTmax,dT);  % est. of max. time used for a parfor iteration
    
    % update others
    idx1 = idxs(end) + 1;
    ntasks_old = ntasks;
end

% close pool
delete(gcp('nocreate'));

end