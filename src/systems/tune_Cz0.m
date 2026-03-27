close all;
opts.sys = 7;
save_flag = false;
rng default;

[sys_dir,~,~] = fileparts(which(mfilename));

[plant,sys_subdir,fn_Cz0] = get_sys_info(opts);
[~,B,C,~,~] = plant2ABCDK(plant);
[ny,nx] = size(C); nu = size(B,2);

switch opts.sys
    case {1,2,3} % for plant from Landau1995
        W1 = makeweight(33,5,0.5);  W1 = c2d(W1,plant.Ts,'tustin');
        W3 = makeweight(0.5,20,20); W3 = c2d(W3,plant.Ts,'tustin');
        W2 = [];
    case 4 % for plant from Bemporad2002
        W1 = makeweight_dB_Hz(80,1,-20,plant.Ts);
        W2 = ss(1e-2);
        W3 = makeweight_dB_Hz(-5,2,20,plant.Ts);
    case 5 % for plant from Favoreel1999
        W1 = makeweight_dB_Hz(20,0.05,-5,plant.Ts);
        W2 = [];
        W3 = makeweight_dB_Hz(-20,0.2,20,plant.Ts);
    case 6 % for plant from Wang2023
        W1 = makeweight_dB_Hz(80,0.05,-5,plant.Ts);
        W2 = [];
        W3 = makeweight_dB_Hz(-30,0.2,5,plant.Ts);
end

%% mixed-sensitivity analysis
G = plant(:,1);
if opts.sys ~= 7
    % if opts.sys = 7 -> use initial controller Cz0 provided by Wang et al. (2023)
    % otherwise, overwrite empty Cz0
    [Cz0,Ms,gamma] = mixsyn(G,W1,W2,W3);
end

switch opts.sys
    
    case 1
        Cz0 = minreal(Cz0);
        nK = 5;
        Cz0 = hinfstruct_wrapper_v2(Cz0,nK,G,plant,W1,W2,W3);
        Cz0 = employ_integrator(Cz0);
        [S,KS,T,Ms] = get_sensitivities(Cz0,G,W1,W2,W3);
        if save_flag
            save(fullfile(sys_dir,'Landau1995','Cz0_Landau1995_D0.mat'),'Cz0');
        end
        
    case 2
        Cz0 = minreal(Cz0);
        nK = 50;
        [Cz0,gamma] = hinfstruct_wrapper_v2(Cz0,nK,G,plant,W1,W2,W3);
        Cz0 = employ_integrator(Cz0);
        [S,KS,T,Ms] = get_sensitivities(Cz0,G,W1,W2,W3);
        if save_flag
            save(fullfile(sys_dir,'Landau1995','Cz0_Landau1995_D0_n50.mat'),'Cz0');
        end
        
    case 3
        Cz0 = minreal(Cz0);
        Cz0 = employ_integrator(Cz0);
        [S,KS,T,Ms] = get_sensitivities(Cz0,G,W1,W2,W3);
        if save_flag
            save(fullfile(sys_dir,'Landau1995','Cz0_Landau1995.mat'),'Cz0');
        end

    case 4
        Cz0 = minreal(Cz0);
        nK = 4;
        [Cz0,gamma] = hinfstruct_wrapper_v2(Cz0,nK,G,plant,W1,W2,W3);
        Cz0 = employ_integrator(Cz0);
        [S,KS,T,Ms] = get_sensitivities(Cz0,G,W1,W2,W3);
        if save_flag
            save(fullfile(sys_dir,'Bemporad2002','Cz0_Bemporad2002.mat'),'Cz0');
        end

    case 5
        Cz0 = minreal(Cz0);
        nK = 7;
        [Cz0,gamma] = hinfstruct_wrapper_v2(Cz0,nK,G,plant,W1,W2,W3);
        % Cz0 = employ_integrator(Cz0); % plant already has an integrator
        [S,KS,T,Ms] = get_sensitivities(Cz0,G,W1,W2,W3);
        if save_flag
            save(fullfile(sys_dir,'Favoreel1999','Cz0_Favoreel1999.mat'),'Cz0');
        end
    
    case 6
        Cz0 = minreal(Cz0);
        nK = 5;
        [Cz0,gamma] = hinfstruct_wrapper_v2(Cz0,nK,G,plant,W1,W2,W3);
        Cz0 = employ_integrator(Cz0);
        [S,KS,T,Ms] = get_sensitivities(Cz0,G,W1,W2,W3);
        if save_flag
            save(fullfile(sys_dir,'Wang2023','Cz0_Wang2023.mat'),'Cz0');
        end

    case 7
        Cz0 = load(fullfile(sys_dir,'Wang2023','Cz0_Wang2023_provided.mat'),'Cz0').Cz0;
        S = feedback(1,G*Cz0);
        KS = Cz0*S;
        T = 1-S;
end

%% closed-loop system
% naming signals
sigs.uk = arrayfun(@(j) sprintf('u0_%d', j), 1:nu, 'UniformOutput', false);
sigs.ek = arrayfun(@(j) sprintf('e0_%d', j), 1:ny, 'UniformOutput', false);
sigs.yk = arrayfun(@(j) sprintf('y0_%d', j), 1:ny, 'UniformOutput', false);
plant.u(1:nu)     =  sigs.uk;
plant.u(nu+1:end) =  sigs.ek;
plant.y = sigs.yk;

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

%% plotting
if opts.sys ~= 7 % skip because Cz0 was loaded, not obtained by tuning
    figure();
    tiledlayout(3,1,'TileSpacing','compact');
    
    ax4 = nexttile;
    if isempty(W2)
        ax41 = sigmaplot(S,'b',KS,'r',T,'g',ss(gamma/W1),'b-.',ss(gamma/W3),'g-.');
        legend('S','KS','T','\gamma/W1','\gamma/W3','Location','SouthWest')
    else
        ax41 = sigmaplot(S,'b',KS,'r',T,'g',ss(gamma/W1),'b-.',ss(gamma/W2),'r-.',ss(gamma/W3),'g-.');
        legend('S','KS','T','\gamma/W1','\gamma/W2','\gamma/W3','Location','SouthWest')
    end
    ax41.FrequencyUnit = 'Hz';
    grid on;
    
    ax5 = nexttile;
    if isempty(W2)
        ax51 = sigmaplot(Ms(1,:),Ms(2,:),Ms);
        Ms_leg = {'$W_1 S_{\infty}$','$W_3 T_{\infty}$','$M_{\infty}$'};
    else
        ax51 = sigmaplot(Ms(1,:),Ms(2,:),Ms(3,:),Ms);
        Ms_leg = {'$W_1 S_{\infty}$','$W_2 KS_{\infty}$','$W_3 T_{\infty}$','$M_{\infty}$'};
    end
    ax51.FrequencyUnit = 'Hz';
    grid on;
    leg = legend(Ms_leg);
    leg.Interpreter = 'latex';
    leg.String = Ms_leg;
    
    ax6 = nexttile;
    bodeopts = bodeoptions;
    bodeopts.PhaseVisible = 'off';
    bodeopts.FreqUnits = 'Hz';
    if isempty(W2)
        ax61 = bodeplot(W1,W3,ss(1),'--',bodeopts);
        legend('W1','W3','0 dB')
    else
        ax61 = bodeplot(W1,W2,W3,ss(1),'--',bodeopts);
        legend('W1','W2','W3','0 dB')
    end
    grid on;
    BodeLines = findall(ax6.Children,'type','line');
    
    linkaxes([ax4 ax5 ax6],'x');
    linkaxes([ax4 ax5],'y');
end

figure;
step(T)

%% Local functions
function W = makeweight_dB_Hz(dc_mag_dB,wc_Hz,hf_mag_dB,Ts)
    dc_mag_abs = db2mag(dc_mag_dB);
    hf_mag_abs = db2mag(hf_mag_dB);
    wc_rad = 2*pi*wc_Hz;
    
    W = makeweight(dc_mag_abs,wc_rad,hf_mag_abs,Ts);
end

% introduce integrator
function Cz0 = employ_integrator(Cz0)
    [tz,po,kg] = zpkdata(Cz0);
    RealPo = real(po{1});
    ImagPo = imag(po{1});
    idxReal = abs(ImagPo) == 0;
    if ~any(idxReal)
        idxReal = abs(ImagPo) < 1e-10;
    end
    if ~any(idxReal)
        fprintf('No integrator introduced because no real poles were found\n');
        return;
    end
    RealPo = RealPo.*idxReal;
    [~,idx_max] = max(RealPo);
    po{1}(idx_max) = 1;            % change largest stable to z=1
    Cz0 = zpk(tz,po,kg,Cz0.Ts);
    Cz0 = ss(tf(Cz0));
end

function [Cz0,gamma] = hinfstruct_wrapper_v1(Cz0,nK,G,plant,W1,W3)
    Cz1 = tunableSS('K',nK,1,1,plant.Ts,'companion');
    % ------- force Cz1.D = 0 -------
    Cz1.D.Value = 0;
    Cz1.D.Free  = false;   % forces Cz1.D = 0
    % ------- initialize rest w/ mixsyn result -----------
    [Cz1.A.Value,Cz1.B.Value,Cz1.C.Value,~] = ssdata(Cz0);

    S = feedback(1,G*Cz1);     % sensitivity
    T = feedback(G*Cz1,1);     % complementary sensitivity
    
    CL = [ W1*S ; W3*T ];
    opt = hinfstructOptions('Display','final');
    [CLopt,gamma] = hinfstruct(CL,opt);

    Cz0 = getBlockValue(CLopt,'K');
end

function [Cz0,gamma] = hinfstruct_wrapper_v2(Cz0,nK,G,plant,W1,W2,W3)
    [num,den] = tfdata(tf(ss(Cz0.A,Cz0.B,Cz0.C,Cz0.D*0,plant.Ts))); num = num{1}; den = den{1};
    nx = size(Cz0.A,1);
    Cz1 = tunableSS('K2',nK,1,1,plant.Ts,'full');
    Cz1.A.Free(1:nK-1,:) = false;
    Cz1.A.Free(end,:) = true;
    Cz1.A.Value(1:nK-1,:) = [zeros(nK-1,1) eye(nK-1)];
    Cz1.B.Free(:) = false; Cz1.B.Value(:,1) = [zeros(nK-1,1); 1];
    Cz1.D.Value = 0; Cz1.D.Free(:) = false;
    % ------- initialize rest w/ mixsyn result -----------
    Cz1.C.Value(1,1:nx) = fliplr(num(1,2:end));
    Cz1.A.Value(end,1:nx) = -fliplr(den(1,2:end));

    S = feedback(1,G*Cz1);     % sensitivity
    T = feedback(G*Cz1,1);     % complementary sensitivity
    KS = Cz1*S;                % control sensitivity
    
    if isempty(W2)
        CL = [ W1*S ; W3*T ];
    else
        CL = [ W1*S ; W2*KS ; W3*T ];
    end
    opt = hinfstructOptions('Display','final');
    [CLopt,gamma] = hinfstruct(CL,opt);

    Cz0 = getBlockValue(CLopt,'K2');
    Cz0 = ss(real(Cz0.A),real(Cz0.B),real(Cz0.C),real(Cz0.D),plant.Ts); % forces real-valued controller
end

function [S,KS,T,Ms] = get_sensitivities(Cz0,G,W1,W2,W3)
    S = feedback(1,G*Cz0);
    KS = Cz0*S;
    T = 1-S;
    if isempty(W2)
        Ms = [W1*S;W3*T];
    else
        Ms = [W1*S;W2*KS;W3*T];
    end
end