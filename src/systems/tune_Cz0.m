close all;
opts.sys = 8;
[plant,sys_subdir,fn_Cz0] = get_sys_info(opts);
[~,B,C,~,~] = plant2ABCDK(plant);
[ny,nx] = size(C); nu = size(B,2);

switch opts.sys
    case {1,5,6,7,8} % for plant from Landau1995
        W1 = makeweight(33,5,0.5);  W1 = c2d(W1,plant.Ts,'tustin');
        W3 = makeweight(0.5,20,20); W3 = c2d(W3,plant.Ts,'tustin');
        W2 = [];
    case 2 % for plant from Bemporad2002
        W1 = makeweight(db2mag(80),[pi/plant.Ts*0.9 1],0.5,plant.Ts);
        W3 = makeweight(0.5,[pi/plant.Ts*0.95 1],20,plant.Ts);
        W2 = ss(1e-1);
    case 3 % for plant from Favoreel1999
        W1 = makeweight(db2mag(80),[pi/plant.Ts*0.58 1],0.85,plant.Ts);
        W3 = makeweight(0.85,[pi/plant.Ts*0.60 1],20,plant.Ts);
        W2 = ss(1e-1);
    case 4 % for plant from Wang2023
        W1 = makeweight(db2mag(80),[pi/plant.Ts*0.58 1],0.85,plant.Ts);
        W3 = makeweight(0.85,[pi/plant.Ts*0.60 1],20,plant.Ts);
        W2 = [];
end

%% mixed-sensitivity analysis
G = plant(:,1);
if opts.sys ~= 9
    % if opts.sys = 9 -> use initial controller Cz0 provided by Wang et al. (2023)
    % otherwise, overwrite empty Cz0
    [Cz0,Ms,gamma] = mixsyn(G,W1,W2,W3);
end

switch opts.sys
    case {1,5}
        Cz0 = minreal(Cz0);
        Cz0 = employ_integrator(Cz0);
        [S,KS,T,Ms] = get_sensitivities(Cz0,G,W1,W3);

    case 6
        Cz0 = minreal(Cz0);
        nK = 5;
        Cz0 = hinfstruct_wrapper_v1(Cz0,nK,G,plant,W1,W3);
        Cz0 = employ_integrator(Cz0);
        [S,KS,T,Ms] = get_sensitivities(Cz0,G,W1,W3);

    case 7
        Cz0 = minreal(Cz0);
        nK = 20;
        Cz0 = hinfstruct_wrapper_v2(Cz0,nK,G,plant,W1,W3);
        [S,KS,T,Ms] = get_sensitivities(Cz0,G,W1,W3);
        
    case 8
        Cz0 = minreal(Cz0);
        nK = 50;
        Cz0 = hinfstruct_wrapper_v2(Cz0,nK,G,plant,W1,W3);
        [S,KS,T,Ms] = get_sensitivities(Cz0,G,W1,W3);

    otherwise
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
figure();
tiledlayout(2,3,'TileSpacing','compact');

ax1 = nexttile;
if isempty(W2)
    ax11 = sigmaplot(S,'b',KS,'r',T,'g',ss(gamma/W1),'b-.',ss(gamma/W3),'g-.');
    legend('S','KS','T','1/W1','1/W3','Location','SouthWest')
else
    ax11 = sigmaplot(S,'b',KS,'r',T,'g',ss(gamma/W1),'b-.',ss(gamma/W2),'r-.',ss(gamma/W3),'g-.');
    legend('S','KS','T','1/W1','1/W2','1/W3','Location','SouthWest')
end
ax11.MagnitudeUnit = 'abs';
ax11.MagnitudeScale = 'linear';
ax11.FrequencyUnit = 'Hz';
grid on;

ax2 = nexttile;
if isempty(W2)
    ax21 = sigmaplot(Ms(1,:),Ms(2,:),Ms);
    Ms_leg = {'$W_1 S_{\infty}$','$W_3 T_{\infty}$','$M_{\infty}$'};
else
    ax21 = sigmaplot(Ms(1,:),Ms(2,:),Ms(3,:),Ms);
    Ms_leg = {'$W_1 S_{\infty}$','$W_2 KS_{\infty}$','$W_3 T_{\infty}$','$M_{\infty}$'};
end
ax21.MagnitudeUnit = 'abs';
ax21.MagnitudeScale = 'linear';
ax21.FrequencyUnit = 'Hz';
leg = legend(Ms_leg);
leg.Interpreter = 'latex';
leg.String = Ms_leg;
grid on;

ax3 = nexttile;

ax4 = nexttile;
if isempty(W2)
    ax41 = sigmaplot(S,'b',KS,'r',T,'g',ss(gamma/W1),'b-.',ss(gamma/W3),'g-.');
else
    ax41 = sigmaplot(S,'b',KS,'r',T,'g',ss(gamma/W1),'b-.',ss(gamma/W2),'r-.',ss(gamma/W3),'g-.');
end
ax41.FrequencyUnit = 'Hz';
grid on;

ax5 = nexttile;
if isempty(W2)
    ax51 = sigmaplot(Ms(1,:),Ms(2,:),Ms);
else
    ax51 = sigmaplot(Ms(1,:),Ms(2,:),Ms(3,:),Ms);
end
ax51.FrequencyUnit = 'Hz';
grid on;

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

axes(ax3);
for k = 1:length(BodeLines)
    Line = BodeLines(k);
    semilogx(Line.XData,db2mag(Line.YData)); hold on;
end
grid on;
xlim(ax3,ax6.XLim);

linkaxes([ax1 ax2 ax3 ax4 ax5 ax6],'x');
linkaxes([ax1 ax2],'y');
linkaxes([ax4 ax5],'y');

figure;
step(T)

%% Helper functions
% introduce integrator
function Cz0 = employ_integrator(Cz0)
    [tz,po,kg] = zpkdata(Cz0);
    [~,idx] = sort(abs(po{1}));
    po{1}(idx(end)) = 1;            % change largest stable to z=1
    Cz0 = zpk(tz,po,kg,Cz0.Ts);
    Cz0 = ss(tf(Cz0));
end

function Cz0 = hinfstruct_wrapper_v1(Cz0,nK,G,plant,W1,W3)
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

function Cz0 = hinfstruct_wrapper_v2(Cz0,nK,G,plant,W1,W3)
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
    
    CL = [ W1*S ; W3*T ];
    opt = hinfstructOptions('Display','final');
    [CLopt,gamma] = hinfstruct(CL,opt);

    Cz0 = getBlockValue(CLopt,'K2');
    Cz0 = ss(real(Cz0.A),real(Cz0.B),real(Cz0.C),real(Cz0.D),plant.Ts); % forces real-valued controller

    % insert an integrator (to complement H-inf design)
    Poles = eig(Cz0.A);
    % Extract only real poles (exclude those with nonzero imaginary part)
    isRealPole = abs(imag(Poles)) == 0 ;%< 1e-10;
    realPoles = Poles(isRealPole);
    [~, idxMax] = max(abs(real(realPoles)));  % Find pole with largest real part
    idxInOriginal = find(isRealPole);
    Poles(idxInOriginal(idxMax)) = 1;    % Set largest real pole to z=1
    denP = poly(diag(Poles));
    [num,den] = tfdata(tf(Cz0));
    Cz0 = ss(tf(num{1},denP,plant.Ts));
end

function [S,KS,T,Ms] = get_sensitivities(Cz0,G,W1,W3)
    S = feedback(1,G*Cz0);
    KS = Cz0*S;
    T = 1-S;
    Ms = [W1*S;W3*T];
end