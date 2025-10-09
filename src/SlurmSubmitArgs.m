function SubmitTxt = SlurmSubmitArgs(jn,Time,opts)
arguments
    jn (1,:) char
    Time (1,1) double {mustBeFinite,mustBeReal,mustBePositive} = 15; % minutes
    opts.part (1,:) char = 'compute';
    opts.accnt (1,:) char = 'research-me-dcsc';
    opts.nodes  = [] %(1,1) double {mustBeFinite,mustBeReal,mustBeInteger,mustBePositive} = 1;                % number of nodes used
    opts.tpn    = [] %(1,1) double {mustBeFinite,mustBeReal,mustBeInteger,mustBeGreaterThan(opts.tpn,1)} = 2; % tasks per node (leave one for submitting script/function)
    opts.ntasks = []
    opts.cpt (1,1)   double {mustBeFinite,mustBeReal,mustBeInteger,mustBePositive} = 1;                % CPUs per task
    opts.GB  (1,1)   double {mustBeFinite,mustBeReal,mustBePositive} = 3.9;                            % GB per CPU

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
else %if ~isempty(opts.tpn)
    check_FinRealPosInt(opts.tpn);
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
    end
end
end

function check_FinRealPosInt(ntasks)
arguments
    ntasks (1,1) double {mustBeFinite,mustBeReal,mustBeInteger,mustBePositive}
end
end