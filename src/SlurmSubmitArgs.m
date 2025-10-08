function SlurmSubmitArgsTxt = SlurmSubmitArgs(jn,Time,opts)
arguments
    jn (1,:) char
    Time (1,1) double {mustBeFinite,mustBeReal,mustBePositive} = 15; % minutes
    opts.part (1,:) char = 'compute';
    opts.nodes (1,1) double {mustBeFinite,mustBeReal,mustBeInteger,mustBePositive} = 1;                % number of nodes used
    opts.tpn (1,1)   double {mustBeFinite,mustBeReal,mustBeInteger,mustBeGreaterThan(opts.tpn,1)} = 2; % tasks per node (leave one for submitting script/function)
    opts.cpt (1,1)   double {mustBeFinite,mustBeReal,mustBeInteger,mustBePositive} = 1;                % CPUs per task
    opts.GB  (1,1)   double {mustBeFinite,mustBeReal,mustBeInteger,mustBePositive} = 3;                % GB per CPU
    opts.accnt (1,:) char = 'research-me-dcsc';
end
Time = char(duration(0,Time,0));
SlurmSubmitArgsTxt = sprintf('--job-name=%s --time=%s --account=%s --partition=%s --nodes=%d --ntasks-per-node=%d --cpus-per-task=%d --mem-per-cpu=%.1fGB',...
                              jn,           Time,     opts.accnt,  opts.part,     opts.nodes, opts.tpn,           opts.cpt,          opts.GB);
end