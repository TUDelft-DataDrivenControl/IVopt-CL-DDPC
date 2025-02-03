function [uf, yf_hat,varargout] = optimizer_solve(obj,opt)
% solve problem using call to yalmip Optimizer object
arguments
    obj
    opt.rf (:,:) double = []
end
if ~isempty(opt.rf)
    obj.rf = opt.rf;
end
switch obj.options.Framework
    case 1 % using YALMIP
        [uf, yf_hat] = yalmip_solver(obj);
    case 2 % using CasADi with Opti
        [uf,yf_hat] = CasADiOpti_solver(obj);
    case 3 % using CasADi without Opti
        [uf,yf_hat,sol_stat,tSolve,PE_stat] = CasADi_solver(obj);
        varargout{1} = sol_stat;
        varargout{2} = tSolve;
        varargout{3} = PE_stat;
end
end

%% Solve with YALMIP optimizer
function [uf, yf_hat] = yalmip_solver(obj)
if obj.options.ExplicitPredictor
    [Lu,Ly,Gu] = obj.getPredictorMatrices();
    [sol,errorcode] = obj.Prob.Optimizer(Lu,Ly,Gu,obj.up,obj.yp,obj.rf);
else
    [sol,errorcode] = obj.Prob.Optimizer(obj.LHS,obj.up,obj.yp,obj.rf);
end
if errorcode~=0 && (any(isnan(sol{1})) || any(isnan(sol{2})))
    error(yalmiperror(errorcode))
end
[uf, yf_hat] = deal(sol{:});
end

%% Solve with CasADi + Opti
function [uf,yf_hat] = CasADiOpti_solver(obj)
if obj.options.ExplicitPredictor
    [Lu,Ly,Gu] = obj.getPredictorMatrices(true);
    [uf, yf_hat] = obj.Prob.Optimizer(Lu,Ly,Gu,obj.up,obj.yp,obj.rf);
else
    [uf,yf_hat] = obj.Prob.Optimizer(obj.LHS,obj.up,obj.yp,obj.rf);
end
uf     = full(uf);
yf_hat = full(yf_hat);
end

%% Solve with CasADi without Opti
function [uf,yf_hat,varargout] = CasADi_solver(obj)
% get parameter vector
if obj.options.ExplicitPredictor
    [obj.Prob.Lu,obj.Prob.Ly,obj.Prob.Gu] = obj.getPredictorMatrices(true); % true to ensure full Gu
    par_vec = obj.Prob.get_p(obj.up,obj.yp,obj.rf,obj.Prob.Lu,obj.Prob.Ly,obj.Prob.Gu);
else
    par_vec = obj.Prob.get_p(obj.up,obj.yp,obj.rf,obj.LHS);
end

% -> initial guess: uf = 0, yf = Lu*up+Ly*yp
yf0 = obj.Prob.Lu*obj.up(:)+obj.Prob.Ly*obj.yp(:);
obj.Prob.res.x  = [zeros(obj.nu*obj.f,1);yf0];

try
    % try to solve regular problem
    obj.Prob.stat = 0;  % <- indicates unsuccessful solve(s)
    solve_timer = tic;
    res = obj.Prob.p2res(par_vec);
    tSolve = toc(solve_timer);
    obj.Prob.res = res; % saving result
    obj.Prob.stat = 1;  % indicating solver used
catch
    try
        % try to solve original problem with softened constraints
        % -> solve relaxed problem
        obj.Prob.stat = 2;
        solve_timer = tic;
        res = obj.Prob.backup.p2res(par_vec);
        tSolve = toc(solve_timer);
        obj.Prob.res  = res;
        obj.Prob.stat = 3;
    catch Error
        disp(['No feasible solution found. Solution status:', num2str(obj.Prob.stat) ])
        error(Error.message);
    end
end
% get uf, yf_hat
[uf, yf_hat] = obj.Prob.res2ufyf(res);

% saving
obj.Prob.yf = yf_hat;
obj.Prob.uf = uf;

% indicate solution status
varargout{1} = obj.Prob.stat;
varargout{2} = tSolve;

% Persistency of excitation - rank check
Z = [obj.Upf; obj.Yp];
if rank(Z) < min(size(Z)) % not full rank
    PE_stat = 0;
else % full rank
    PE_stat = 1;
end
varargout{3} = PE_stat;

end