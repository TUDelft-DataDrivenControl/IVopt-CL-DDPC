function [plant,nu,ny,Cz0,Tcl0,opts,sigs] = init_sims(opts)
%% initializes closed-loop simulations by getting necessary variables

%% choose a  plant model
switch opts.sys
    case {1,5,6,7,8}
        [plant,nx,nu,ny,A,B,C,D,K,~] = model_Landau1995();
        W1 = makeweight(33,5,0.5);  W1 = c2d(W1,plant.Ts,'tustin');
        W3 = makeweight(0.5,20,20); W3 = c2d(W3,plant.Ts,'tustin');
        
        switch opts.sys
            case 1
                fn_Cz0 = 'Cz0_Landau1995.mat';
            case 5
                fn_Cz0 = 'Cz0_Landau1995.mat'; % still tuned for opts.sys = 1
                % change system such that it has no input-output delay
                [b,a] = ss2tf(A,B(:,1),C,D(:,1));
                b2 = circshift(b,-3);
                plant = minreal([tf(b2,a,plant.Ts) plant(:,2)]);
                A = plant.A;
                B = plant.B(:,1);
                C = plant.C;
                D = plant.D(:,1);
            case 6
                fn_Cz0 = 'Cz0_Landau1995_D0.mat';
            case 7
                fn_Cz0 = 'Cz0_Landau1995_D0_n20.mat';
            case 8
                fn_Cz0 = 'Cz0_Landau1995_D0_n50.mat';
        end
    case 2
        [plant,nx,nu,ny,A,B,C,D,K,~] = model_Bemporad2002(At_poles=[0.95, 0.9]);
        plant.Ts = 1;
        W1 = makeweight(db2mag(80),[pi/plant.Ts*0.9 1],0.5,plant.Ts);
        W3 = makeweight(0.5,[pi/plant.Ts*0.95 1],20,plant.Ts);
        fn_Cz0 = 'Cz0_Bemporad2002.mat';
    case 3
        % system from Favoreel 1999; SPC: Subspace Predictive Control
        A = [ 4.40 1 0 0 0;
             -8.09 0 1 0 0;
              7.83 0 0 1 0;
             -4.00 0 0 0 1;
              0.86 0 0 0 0];
        B = [0.00098 0.01299 0.01859 0.0033 -0.00002].';
        C = eye(1,5);
        D = 0;
        K = [2.3 -6.64 7.515 -4.0146 0.86336].';
        
        nx = size(A,1); nu = size(B,2); ny = size(C,1);
        
        plant = ss(A,[B K], C, [D eye(ny,nu)],1);
        W1 = makeweight(db2mag(80),[pi/plant.Ts*0.58 1],0.85,plant.Ts);
        W3 = makeweight(0.85,[pi/plant.Ts*0.60 1],20,plant.Ts);
        fn_Cz0 = 'Cz0_Favoreel1999.mat';
    case 4
        [plant,Cz0,nx,nu,ny,A,B,C,D,K,Re] = model_Wang2023();
        fn_Cz0 = 'Cz0_Wang2023.mat';
        plant.Ts = 1;
        W1 = makeweight(db2mag(80),[pi/plant.Ts*0.58 1],0.85,plant.Ts);
        W3 = makeweight(0.85,[pi/plant.Ts*0.60 1],20,plant.Ts);
        W2 = [];%ss(1e-1);
end

% naming signals
sigs.uk = arrayfun(@(j) sprintf('u0_%d', j), 1:nu, 'UniformOutput', false);
sigs.ek = arrayfun(@(j) sprintf('e0_%d', j), 1:ny, 'UniformOutput', false);
sigs.yk = arrayfun(@(j) sprintf('y0_%d', j), 1:ny, 'UniformOutput', false);
plant.u(1:nu)     =  sigs.uk;
plant.u(nu+1:end) =  sigs.ek;
plant.y = sigs.yk;

% saving to opts structure
[opts.ny,opts.nu,opts.nx] = deal(ny,nu,nx);

%% =============== for initial closed-loop simulation ======================
% ----------------- make/load initial controller (Cz0) --------------------
if isfile(fn_Cz0)
    Cz0 = load(fn_Cz0).Cz0;
else
    [Cz0,~,~,~] = mixsyn(plant(:,1:nu),W1,[],W3);
    save(fn_Cz0,"Cz0");
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