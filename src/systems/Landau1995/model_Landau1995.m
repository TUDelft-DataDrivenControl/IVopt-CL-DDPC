function plant = model_Landau1995()
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

%% system matrices
[A,B,C,D] = ssdata(ss(...
            tf([0 0 0 0.28261 0.50666],[1 -1.41833 1.58939 -1.31608 0.88642],...
               [],'Variable', 'z^-1')));
K = 1.9*[0.0939; -0.3433; 0.1063; 1.2058]; % <- not in [1], so from [2]
% Re = 4.81*1e-3 % <- from [2], but unused since it will be set by user

nu = size(B,2);
ny = size(C,1);

%% create plant model
plant = ss(A,[B K], C, [D eye(ny,nu)],0.05);
end