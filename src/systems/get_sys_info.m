function [plant,sys_subdir,fn_Cz0] = get_sys_info(opts)
%% this function gets plant and initial controller information for the system specified in opts.sys
% plant: state-space model of the plant of the form
%           x_{k+1} = A x_k + [B K] [u_k; e_k]      i.e. plant.B = [B K]
%               y_k = C x_k + [D I] [u_k; e_k]      i.e. plant.D = [D I]
% sys_subdir:   subdirectory name in src/systems/ where system and initial controller
%               files are located
% fn_Cz0:       filename of the initial controller (Cz0)

switch opts.sys
    case {1,6,8}
        plant = model_Landau1995();
        sys_subdir = 'Landau1995';
        switch opts.sys
            case 1
                % controller with direct feedthrough, 5 states
                % -> CL-SPC and Transient Predictor now employ biased estimates, added value of an IV all the more apprarent
                fn_Cz0 = 'Cz0_Landau1995.mat';
            case 6
                % controller without direct feedthrough, 5 states
                % -> case shown in paper
                fn_Cz0 = 'Cz0_Landau1995_D0.mat';
            case 8
                % controller without direct feedthrough, 50 states
                % -> interesting to compare performance of IV4 & IV5 (with vs. without Cz0 knowledge)
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