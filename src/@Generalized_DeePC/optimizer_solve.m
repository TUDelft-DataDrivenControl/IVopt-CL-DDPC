function [uf, yf_hat,varargout] = optimizer_solve(obj,opt)
% solve problem using call to optimizer
arguments
    obj
    opt.rf (:,:) double = []
end
if ~isempty(opt.rf)
    obj.rf = opt.rf;
end

if obj.options.ExplicitPredictor
    % get parameter vector
    [obj.Prob.Lu,obj.Prob.Ly,obj.Prob.Gu] = obj.getPredictorMatrices(true); % true to ensure full Gu
    par_vec = obj.Prob.get_p(obj.up,obj.yp,obj.rf,obj.Prob.Lu,obj.Prob.Ly,obj.Prob.Gu);

    % -> initial guess: uf = 0, yf = Lu*up+Ly*yp
    yf0 = obj.Prob.Lu*obj.up(:)+obj.Prob.Ly*obj.yp(:);
    obj.Prob.res.x  = [zeros(obj.nu*obj.f,1);yf0];
else
    % get parameter vector
    par_vec = obj.Prob.get_p(obj.up,obj.yp,obj.rf,obj.LHS);
    
    % -> initial guess: uf = 0, Adyn*[uf(:);yf(:);G(:)] = ulba_dyn
    Adyn = full(obj.Prob.get_Adyn(par_vec));
    Adyn = Adyn(:,obj.nu*obj.f+1:end);
    ulba_dyn = full(obj.Prob.get_uba(par_vec));
    ulba_dyn = ulba_dyn(end-size(Adyn,1)+1:end,1);
    sol1 = pinv(Adyn)*ulba_dyn;
    obj.Prob.res.x = [zeros(obj.nu*obj.f,1);sol1];
end

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
% indicate solution status
varargout{1} = obj.Prob.stat;
varargout{2} = tSolve;

% get uf, yf_hat
[uf, yf_hat] = obj.Prob.res2ufyf(res);

% saving
obj.Prob.yf = yf_hat;
obj.Prob.uf = uf;

% Persistency of excitation - rank check
Z = [obj.Upf; obj.Yp];
if rank(Z) < min(size(Z)) % not full rank
    PE_stat = 0;
else % full rank
    PE_stat = 1;
end
varargout{3} = PE_stat;
end