close all;
opts.sys = 4;
switch opts.sys
    case 1
        [plant,nx,nu,ny,A,B,C,D,K,~] = model_Landau1995();
        ref_fac = 1;
        W1 = makeweight(33,5,0.5);  W1 = c2d(W1,plant.Ts,'tustin');
        W3 = makeweight(0.5,20,20); W3 = c2d(W3,plant.Ts,'tustin');
        W2 = [];
    case 2
        [plant,nx,nu,ny,A,B,C,D,K,~] = model_Bemporad2002(At_poles=[0.95, 0.9]);
        ref_fac = 1;
        plant.Ts = 1;
        W1 = makeweight(db2mag(80),[pi/plant.Ts*0.9 1],0.5,plant.Ts);
        W3 = makeweight(0.5,[pi/plant.Ts*0.95 1],20,plant.Ts);
        W2 = ss(1e-1);
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
        ref_fac = 1e0;
        W1 = makeweight(db2mag(80),[pi/plant.Ts*0.58 1],0.85,plant.Ts);
        W3 = makeweight(0.85,[pi/plant.Ts*0.60 1],20,plant.Ts);
        W2 = ss(1e-1);
    case 4
        [plant,Cz0,nx,nu,ny,A,B,C,D,K,Re] = model_Wang2023();
        plant.Ts = 1;
        W1 = makeweight(db2mag(80),[pi/plant.Ts*0.58 1],0.85,plant.Ts);
        W3 = makeweight(0.85,[pi/plant.Ts*0.60 1],20,plant.Ts);
        W2 = [];%ss(1e-1);
end
%% mixed-sensitivity analysis
G = plant(:,1);
[Cz0,Ms,gamma] = mixsyn(G,W1,W2,W3); gamma
S = feedback(1,G*Cz0);
KS = Cz0*S;
T = 1-S;

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