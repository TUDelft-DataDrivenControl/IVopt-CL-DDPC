function SubmitTxt = SlurmSubmitArgs(jn,Time,opts)
arguments
    jn (1,:) char % job name
    Time (1,1) double {mustBeFinite,mustBeReal,mustBePositive} = 15; % max. runtime in minutes
    opts.part char  = ''; % partition
    opts.accnt char = ''; % account
    opts.nodes  = [] % number of nodes used
    opts.tpn    = [] % tasks per node
    opts.ntasks = [] % number of tasks (single node)
    opts.cpt (1,1)   double {mustBeFinite,mustBeReal,mustBeInteger,mustBePositive} = 1; % CPUs per task
    opts.GB  (1,1)   double {mustBeFinite,mustBeReal,mustBePositive} = 3.9;             % GB per CPU
    opts.out char = ''; % output file for stdout
    opts.err char = ''; % error file for stderr
end
%% futher parsing
if (isempty(opts.tpn) && isempty(opts.ntasks)) || ~isempty(opts.tpn) && ~isempty(opts.ntasks)
    error(sprintf('Specify either, but not both of the below:\n\t1) number of tasks per node: \t     tpn= ...\n\t2) number of tasks (single node): ntasks= ...'));
end
if isempty(opts.ntasks) && isempty(opts.nodes)
    opts.nodes = 1;
end
if ~isempty(opts.ntasks)
    check_FinRealPosInt(opts.ntasks);
else
    check_FinRealPosInt(opts.tpn);
end
if isempty(opts.part) || isempty(opts.accnt)
    % load default settings saved to SlurmSettings.mat
    [src_util_dir, ~  , ~] = fileparts(which(mfilename)); % find src/util directory
    SlurmSettings = load(fullfile(src_util_dir,'..','SlurmSettings.mat'));
    if isempty(opts.part)
        opts.part = SlurmSettings.partition;
    end
    if isempty(opts.accnt)
        opts.accnt = SlurmSettings.account;
    end
end

% remove empty fields
fns = fieldnames(opts);
fns_empty = fns(structfun(@isempty,opts));
opts = rmfield(opts,fns_empty);

% format time
Time = char(duration(0,Time,0));

%% create submission string
SubmitTxt = sprintf('--job-name=%s --partition=%s --time=%s', jn, opts.part, Time);

fns = fieldnames(opts);
for k= 1:length(fns)
    fn = fns{k};
    switch fn
        case 'nodes'
            SubmitTxt = sprintf('%s --nodes=%d', SubmitTxt, opts.nodes);
        case 'ntasks'
            SubmitTxt = sprintf('%s --ntasks=%d', SubmitTxt, opts.ntasks);
        case 'tpn'
            SubmitTxt = sprintf('%s --ntasks-per-node=%d', SubmitTxt, opts.tpn);
        case 'cpt'
            SubmitTxt = sprintf('%s --cpus-per-task=%d', SubmitTxt, opts.cpt);
        case 'GB'
            SubmitTxt = sprintf('%s --mem-per-cpu=%dMB', SubmitTxt, floor(opts.GB*1e3));
        case 'out'
            SubmitTxt = sprintf('%s --output=%s', SubmitTxt, opts.out);
        case 'err'
            SubmitTxt = sprintf('%s --error=%s', SubmitTxt, opts.err);
    end
end
end

function check_FinRealPosInt(ntasks)
arguments
    ntasks (1,1) double {mustBeFinite,mustBeReal,mustBeInteger,mustBePositive}
end
end