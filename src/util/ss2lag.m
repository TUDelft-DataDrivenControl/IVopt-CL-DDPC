function lag = ss2lag(Cz0)
[Ac,~,Cc,~] = ssdata(Cz0);
nu = size(Cc,1);
nxc = size(Ac,1);
Gamma_c = make_ext_obsv(Ac,Cc,nxc);
if rank(Gamma_c)~= nxc
    error('Controller is not observable');
end
% rho = lag of the controller
lag = nan;
for k = 1:nxc
    if rank(Gamma_c(1:nu*k,:)) == nxc
        lag = k;
        break;
    end
end
end