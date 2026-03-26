function [plant,nx,nu,ny,A,B,C,D,K,Cz0] = model_Wang2023(options)
%% System from Wang et al. (2023):
%
% [1] Y. Wang, Y. Qiu, M. Sader, D. Huang, and C. Shang, “Data-Driven
%     Predictive Control Using Closed-Loop Data: An Instrumental Variable
%     Approach,” IEEE Control Systems Letters, vol. 7, pp. 3639–3644, 2023,
%     doi: 10.1109/LCSYS.2023.3340444.
%
% The above authors of [1] also provide an initial controller Cz0, which is also provided below.
% Since K is not specified, it is calculated by pole placement of A-KC, with poles specified by user

arguments
   options.K (3,1) double  = nan(3,1)               % <- not mentioned in [1]
   options.At_poles (1,3) double = [0.95;0.89;0.88] % <- default to calculate K with
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
    K = place(A.',C.',options.At_poles).';
end

nx = 3;
nu = 1;
ny = 1;

%% controller <- acts on error: e_k = r_k - y_k
Ac = [ 1.0257  0.1540;...
      -0.0043  0.9743];
Bc = [-0.2046; -0.0963];
Cc = [-0.4604 -1.5539];
Dc = -0.03;

%% create systems
plant = ss(A,[B K], C, [D eye(ny,nu)],[]);
Cz0 = ss(Ac,Bc,Cc,Dc,[]);
end
