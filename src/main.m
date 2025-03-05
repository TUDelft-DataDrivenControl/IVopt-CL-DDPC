%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           DDPC using an Optimal-IV
%           Authors: R. Dinkla, T. Oomen, J.W. van Wingerden
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;
close all;
clc;
rng default;

%% controller settings

%% user defined constraints
y_max = 1000;
u_max = 15;
du_max = 3.75;

%% make used stochastic data
% -> save now s.t. if there is an error this data is available for debugging
OL_sim_steps = num_steps-CL_sim_steps;
e     = mvnrnd(zeros(ny,1),Re, num_steps).';    % noise realization
u_OL  = mvnrnd(zeros(nu,1),Ru, OL_sim_steps).'; % OL input
du_CL = mvnrnd(zeros(nu,1),Rdu,CL_sim_steps).'; % CL input

% splitting e into OL & CL parts
e_OL = e(:,1:OL_sim_steps);
e_CL = e(:,OL_sim_steps+1:end);

%% initial open loop simulation

% simulate system
[y_OL,~,x_OL] = lsim(plant,[u_OL;e_OL],[],x0);
y_OL = y_OL.'; x_OL = x_OL.';
x0_CL = plant.A*x_OL(:,end) + plant.B*[u_OL(:,end); e_OL(:,end)];

% select used OL data
u_ol = u_OL(:,end-Nbar+1:end);
y_ol = y_OL(:,end-Nbar+1:end);

%% closed-loop operation

% step-simulate system
simulate_step = @(u_k,e_k,x_k) deal(plant.C*x_k + plant.D*[u_k;e_k],... [y_k,
                                    plant.A*x_k + plant.B*[u_k;e_k]);  % x_next]
stage_cost    = @(u,y_k) y_k.'*Q*y_k + u(:,2).'*R*u(:,2) + (u(:,2)-u(:,1)).'*dR*(u(:,2)-u(:,1));

% user defined constraints
con = struct();
con.y_max  =  y_max;
con.y_min  = -y_max;
con.u_max  =  u_max;
con.u_min  = -u_max;
con.du_max = du_max;

% initialize data arrays
Cz = cell(3,1);
[u_CL,y_CL,x_CL,Cost,stat,tSolves] = deal(Cz);
[eLu,eLy,eGu,eObX,PE_stat] = deal(cell(2,1));

% 1) DeePC with IV
Cz{1} =    DeePC(u_ol,y_ol,p,f,N_OL,Q,R,dR,constr=con);

% 2) CL-DeePC with IV
Cz{2} = CL_DeePC(u_ol,y_ol,p,f,N_CL,Q,R,dR,constr=con,EstimateD=false);

for k_c = 1:2
    % initialize data for run with controller
    x_CLr = nan(nx,CL_sim_steps+1); x_CLr(:,1) = x0_CL;
    u_CLr = nan(nu,CL_sim_steps);
    y_CLr = nan(ny,CL_sim_steps);
    cost  = nan(1,CL_sim_steps);
    stat_r= cost;
    tSolves_r = cost;
    PE_stat_r = cost;
    [eLu_r,eLy_r,eGu_r,eObX_r] = deal(cost);

    % get first CL input
    [uf_k,~,stat_r(1),tSolves_r(1), PE_stat_r(1)] = Cz{k_c}.solve(rf=ref(:,1:f));

    % add disturbance input
    u_CLr(:,1) = uf_k(:,1)+du_CL(:,1);
    % simulate step
    [y_CLr(:,1),x_CLr(:,2)] = simulate_step(u_CLr(:,1),e_CL(:,1),x_CLr(:,1));
    
    % analysis
    cost(1) = stage_cost([u_OL(:,end) u_CLr(:,1)],y_CLr(:,1)); % stage cost

    for k = 2:CL_sim_steps
        % get input
        try
            if k_c < 3
                [uf_k,~,stat_r(k),tSolves_r(k), PE_stat_r(k)] = Cz{k_c}.step( u_CLr(:,k-1), y_CLr(:,k-1), rf=ref(:,k:k+f-1));
            else
                [uf_k,~,stat_r(k),tSolves_r(k)] = Cz{k_c}.solve(x_CLr(:,k),   u_CLr(:,k-1), rf=ref(:,k:k+f-1));
            end
        catch Error
            disp(['k_var =',num2str(k_var),'; k_e = ',num2str(k_e),' k1 = ',num2str(k)]);
            error(Error.message)
        end
        % add disturbance input
        u_CLr(:,k) = uf_k(:,1)+du_CL(:,k);
        % simulate step
        [y_CLr(:,k),x_CLr(:,k+1)] = simulate_step(u_CLr(:,k),e_CL(:,k),x_CLr(:,k));
        % calculate cost
        cost(k) = stage_cost(u_CLr(:,k-1:k),y_CLr(:,k));      % stage cost
    end

    % saving CL data from run with controller
    u_CL{k_c} = u_CLr;
    y_CL{k_c} = y_CLr;
    x_CL{k_c} = x_CLr;
    Cost{k_c} = cost;
    if k_c < 3
        eLu{k_c}  = eLu_r;
        eLy{k_c}  = eLy_r;
        eGu{k_c}  = eGu_r;
        eObX{k_c} = eObX_r;
        PE_stat{k_c} = PE_stat_r;
    end
    stat{k_c} = stat_r;
    tSolves{k_c} = tSolves_r;

end % end for k_c
system("echo Saving loop var results");
% save(save_str,'y_OL','x_OL','u_CL','y_CL','x_CL','Cost','eLu','eLy','eGu','eObX','