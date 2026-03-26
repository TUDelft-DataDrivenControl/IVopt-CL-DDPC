function [plant,nx,nu,ny,A,B,C,D,K] = model_Bemporad2002(options)
%% System from Bemporad et al. (2002):
% Nominal plant model from [1], K is determined by pole placement of A-KC or by 
% solving a discrete-time algebraic Riccati equation based user-specified noise properties
%
% [1] A. Bemporad, M. Morari, V. Dua, and E. N. Pistikopoulos, “The explicit
%     linear quadratic regulator for constrained systems,” Automatica, vol. 38,
%     no. 1, pp. 3–20, Jan. 2002, doi: 10.1016/S0005-1098(01)00174-1.

arguments
   options.At_poles (1,2) double = nan(1,2);
   options.Qw       (2,2) double = nan(2,2);
   options.Rv       (1,1) double = nan(1,1);
   options.Swv      (2,1) double = nan(2,1);
end
[Poles,Qw,Rv,Swv] = deal(options.At_poles,options.Qw,options.Rv,options.Swv);

%% system matrices
A = [0.7326 -0.0861; 0.1722 0.9909];
B = [0.0609; 0.0064];
C = [0 1.4142];
D = 0;
K = zeros(2,1);

nx = size(A,1); % nx = 2
nu = size(B,2); % nu = 1
ny = size(C,1); % ny = 1

if all(~isnan(Poles))
    % since nx = 2, calculate gain K analytically:
    p1 = Poles(1); p2 = Poles(2);
    a = A(1,1); b = A(1,2); c = A(2,1); d = A(2,2);
    g = C(1,2);
    k2 = (a-p1-p2+d)/g; d2 = d - k2*g;
    b2 = ( (p1-p2)^2 - (a+d2)^2 )/(4*c) + a*d2/c;
    k1 = (b - b2)/g;
    K = [k1; k2];
elseif all(~isnan(Qw)) && all(~isnan(Rv))
    if any(isnan(Swv))
        [P,K,~] = idare(A.',B.',Qwv,Rv,[],[]);
    else
        [P,K,~] = idare(A.',B.',Qwv,Rv,Swv,[]);
    end
    K = K.';
    % Re = C*P*C.'+Rv; % <- not used since it will be set by user
end

%% create plant model
plant = ss(A,[B K], C, [D eye(ny,nu)],[]);
end