function [plant,nx,nu,ny,A,B,C,D,K,Re] = model_Landau1995(options)
% nominal model from [1], noise properties from [2]
% [1] I. D. Landau, D. Rey, A. Karimi, A. Voda, and A. Franco, “A Flexible
%     Transmission System as a Benchmark for Robust Digital Control,”
%     European Journal of Control, vol. 1, no. 2, pp. 77–96, Jan. 1995,
%     doi: 10.1016/S0947-3580(95)70011-5.
% 
% [2] A. Chiuso, M. Fabris, V. Breschi, and S. Formentin, “Harnessing
%     uncertainty for a separation principle in direct data-driven
%     predictive control,” Automatica, vol. 173, p. 112070, Mar. 2025,
%     doi: 10.1016/j.automatica.2024.112070.

arguments
   options.K (4,1) double  = [0.1784; -0.6523; 0.2020; 2.2910] % <- not in [1], so from [2]
   options.Re (1,1) double = 4.81*1e-3                         % <- varies in [1], used [2]
   options.At_poles (1,4) double = nan(1,4)                    % <- not in [1]
end

%% system matrices
[A,B,C,D] = ssdata(ss(...
            tf([0 0 0 0.28261 0.50666],[1 -1.41833 1.58939 -1.31608 0.88642],...
               [],'Variable', 'z^-1')));
K = options.K;
Re = options.Re;
if all(~isnan(options.At_poles))
    K = place(A.',C.',options.At_poles).';
end

nx = size(A,1);
nu = size(B,2);
ny = size(C,1);

%% create plant models

plant = ss(A,[B K], C, [D eye(ny,nu)],0.05);
end