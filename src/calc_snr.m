function SNR = calc_snr(plant,opts)
%% Calculates the signal-to-noise ratio (SNR)
% Input arguments:
% - plant: ss(A,[B K],C,[D I],[])
% - optional arguments: either (1) or (2)
%   (1) inputs u and innovation noise e
%   (2) input variance Ru and innovation noise variance Re
%
% Output: SNR [-]    <-- not in dB
arguments
    plant           ss
    opts.Ru         double {mustBeMatrix}
    opts.Re         double {mustBeMatrix}
    opts.Nsim (1,1) double
    opts.u          double {mustBeMatrix}
    opts.e          double {mustBeMatrix}
end
%% Parsing arguments
if isfield(opts,'Ru') && isfield(opts,'Re') && isfield(opts,'Nsim')
    [Ru,Re,Nsim] = deal(opts.Ru,opts.Re,opts.Nsim);
    [nu,nu2] = size(Ru);
    [ny,ny2] = size(Re);
    if nu~=nu2 || ny~=ny2
        error('Ru and Ry must be square');
    end
    u = mvnrnd(zeros(nu,1),Ru,Nsim).';
    e = mvnrnd(zeros(ny,1),Re,Nsim).';
elseif isfield(opts,'u') && isfield(opts,'e')
    [u,e] = deal(opts.u,opts.e);
    Nsim  = size(u,2);
    Nsim2 = size(e,2);
    if Nsim ~= Nsim2
        error('Inputs u (nu|Nsim) and noise e (ny|Nsim) must have the same number of columns');
    end
else
    error(['Provide either inputs u (nu|Nsim) and outputs (ny|Nsim) or',...
        'input variance Ru (nu|nu) and innovation noise variance Re (ny|ny)']);
end

%% Simulating system driven by either inputs or noise
y_e0 = lsim(plant,[u;0*e],[]).'; % driven by input: e_k = 0
y_u0 = lsim(plant,[0*u;e],[]).'; % driven by noise: u_k = 0

%% Calculate SNR
varYe0 = SampleVar(y_e0); % driven by input
varYu0 = SampleVar(y_u0); % driven by noise
SNR = 1+varYe0/varYu0;

%% Helper functions
function VarX = SampleVar(x_all)
    [n,N] = size(x_all);
    VarX = var(x_all,1,2);
    if n ~= 1
        VarX = VarX* (N-1)/(N-n); % correction for multi-dimensionality
    end
end

end