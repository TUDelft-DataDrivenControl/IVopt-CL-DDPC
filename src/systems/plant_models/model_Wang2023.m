function [plant,Cz,nx,nu,ny,A,B,C,D,K,Re] = model_Wang2023(options)
% models from
% Wang, Y., Qiu, Y., Sader, M., Huang, D., & Shang, C. (2023).
% Data-Driven Predictive Control Using Closed-Loop Data: An Instrumental Variable Approach.
% IEEE Control Systems Letters, 7, 3639–3644.

arguments
   options.K (3,1) double  = nan(3,1)   % <- not specified!
   options.Re (1,1) double = 1          % <- varies in paper
   options.At_poles (1,3) double = nan(1,3) % <- not mentioned in paper
end
%% plant model
A = [0.9261  0.0534  0.0382; ...
     0.0907  0.9322 -0.0161; ...
     0.0939  0.0139  0.8753];
B = [0.3609; 0.3064; -0.0918];
C = [-0.6995 0.91    -0.1381];
D = 1;
K = options.K;
if any(isnan(K))
    if any(isnan(options.At_poles))
        options.At_poles = [0.95;0.89;0.88];
    end
    K = place(A.',C.',options.At_poles).';
end
Re = options.Re;

nx = 3;
nu = 1;
ny = 1;

%% controller <- acts on error: e_k = r_k-y_k
Ac = [ 1.0257  0.1540;...
      -0.0043  0.9743];
Bc = [-0.2046; -0.0963];
Cc = [-0.4604 -1.5539];
Dc = -0.03;

%% create systems
plant = ss(A,[B K], C, [D eye(ny,nu)],[]);
Cz = ss(Ac,Bc,Cc,Dc,[]);
end
