% model from
% Bemporad, Alberto, Morari, Manfred, Dua, Vivek, & Pistikopoulos, Efstratios N. (2002).
% The explicit linear quadratic regulator for constrained systems. Automatica, 38(1), 3–20.
A = [0.7326 -0.0861; 0.1722 0.9909];
B = [0.0609; 0.0064];
C = [0 1.4142];
D = 0;
K = zeros(2,1);
Re = 0;

nx = size(A,1);
nu = size(B,2);
ny = size(C,1);

%% create plant models

plant = ss(A,[B K], C, [D eye(ny,nu)],[]);