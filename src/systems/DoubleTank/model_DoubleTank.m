% model from page 153 of the book "Filtering and System Identification" by
% Michel Verhaegen and Vincent Verdult (2007)
A = [0.9512 0; 0.0476 0.9512];
B = [0.0975; 0.0024];
C = [0 1];
D = 0;
Bw = [0.0975 0; 0.0024 0.0975];
Dv = 1;

Rv = 0.0125;
Svw= [0;0.005];
Qw = 0.01*eye(2);

% stationary innovation noise properties
[P,K,L] = idare(A.',C.',Qw,Rv,Svw,[]);
K = K.';
Re = C*P*C.'+Rv;

nx = size(A,1);
nu = size(B,2);
ny = size(C,1);

%% create plant models

% Kalman filter
plant = ss(A,[B K], C, [D eye(ny,nu)],[]);