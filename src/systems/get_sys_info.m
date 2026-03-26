function [plant,sys_subdir,fn_Cz0] = get_sys_info(opts)
%% this function gets plant and initial controller information for the system specified in opts.sys
% plant: state-space model of the plant of the form
%           x_{k+1} = A x_k + [B K] [u_k; e_k]      i.e. plant.B = [B K]
%               y_k = C x_k + [D I] [u_k; e_k]      i.e. plant.D = [D I]
% sys_subdir:   subdirectory name in src/systems/ where system and initial controller
%               files are located
% fn_Cz0:       filename of the initial controller (Cz0)

switch opts.sys
    case {1,5,6,7,8}
        plant = model_Landau1995();
        sys_subdir = 'Landau1995';
        switch opts.sys
            case 1
                fn_Cz0 = 'Cz0_Landau1995.mat';
            case 5
                fn_Cz0 = 'Cz0_Landau1995.mat';
                % modify plant such that it has no input-output delay
                [A,B,C,D,~] = plant2ABCDK(plant);
                [b,a] = ss2tf(A,B(:,1),C,D(:,1));
                b2 = circshift(b,-3);
                plant = minreal([tf(b2,a,plant.Ts) plant(:,2)]);
            case 6
                fn_Cz0 = 'Cz0_Landau1995_D0.mat';
            case 7
                fn_Cz0 = 'Cz0_Landau1995_D0_n20.mat';
            case 8
                fn_Cz0 = 'Cz0_Landau1995_D0_n50.mat';
        end
    case 2
        plant = model_Bemporad2002(At_poles=[0.95, 0.9]);    
        sys_subdir = 'Bemporad2002';
        fn_Cz0 = 'Cz0_Bemporad2002.mat';
    case 3
        plant = model_Favoreel1999();
        sys_subdir = 'Favoreel1999';
        fn_Cz0 = 'Cz0_Favoreel1999.mat';
    case {4,9}
        plant = model_Wang2023();
        sys_subdir = 'Wang2023';
        if opts.sys == 4
            fn_Cz0 = 'Cz0_Wang2023.mat';
        else % opts.sys = 9
            fn_Cz0 = 'Cz0_Wang2023_provided.mat';
        end
end