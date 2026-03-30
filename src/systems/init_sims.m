function [plant,nu,ny,Cz0,Tcl0,opts,sigs] = init_sims(opts)
%% function that gets variables with which to initialize simulations

%% get system information: plant & initial controller
[plant,sys_subdir,fn_Cz0] = get_sys_info(opts);
[~,B,C,~,~] = plant2ABCDK(plant);
[ny,nx] = size(C); nu = size(B,2);

% naming signals
sigs.uk = arrayfun(@(j) sprintf('u0_%d', j), 1:nu, 'UniformOutput', false);
sigs.ek = arrayfun(@(j) sprintf('e0_%d', j), 1:ny, 'UniformOutput', false);
sigs.yk = arrayfun(@(j) sprintf('y0_%d', j), 1:ny, 'UniformOutput', false);
plant.u(1:nu)     =  sigs.uk;
plant.u(nu+1:end) =  sigs.ek;
plant.y = sigs.yk;

% saving to opts structure
[opts.ny,opts.nu,opts.nx] = deal(ny,nu,nx);

%% =============== for initial closed-loop simulation =====================
% ----------------- load initial controller (Cz0) -------------------------
fn_Cz0 = fullfile('systems',sys_subdir,fn_Cz0);
if isfile(fn_Cz0)
    % overwrite empty initial controller Cz0
    Cz0 = load(fn_Cz0).Cz0;
else
    error('Initial controller Cz0 not found: %s\nRun the corresponding tune_Cz0_*.m script to generate it first.', fn_Cz0);
end

% naming signals
Cz0.u = arrayfun(@(j) sprintf('er0_%d', j), 1:ny, 'UniformOutput', false);
Cz0.y = plant.u(1:nu);

% ----------------- make initial closed-loop system (Tcl0) ----------------
fbsum = cell(ny,1);
for k = 1:ny
    fbsum{k} = sumblk(sprintf('er0_%d = r0_%d - y0_%d', k,k,k));
end
sigs.rk = arrayfun(@(j) sprintf('r0_%d', j), 1:ny, 'UniformOutput', false); % r_k
conOpts = connectOptions("Simplify",false);
Tcl0 = connect(Cz0,plant,fbsum{:},[sigs.rk(:).',sigs.ek(:).'],[sigs.uk(:).',sigs.yk(:).'],conOpts);

end