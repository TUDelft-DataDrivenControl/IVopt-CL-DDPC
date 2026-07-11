function run_Re(ii,opts2,spRe,nRe,seeds,Re_all,Ncl,Nbar,ny,plant,subdir1,sigs,Cz0,Tcl0,yr0,yr1,ur1,proj_dir)
%% This function acts as a wrapper for main_MC.m for a single value of Re

opts = opts2;
s = struct;

% ------------------------------ define Re run ----------------------------
[~,iRe] = ind2sub([spRe,nRe],ii);
Re = Re_all(iRe);
opts.Re = Re;

% ------------------------- define seed run & noise -----------------------
seed = seeds(ii); opts.seed = seed;
rng(seed);
s.e0 = mvnrnd(zeros(ny,1),Re,Nbar).'; % innovation noise
s.e1 = mvnrnd(zeros(ny,1),Re,Ncl).';  % innovation noise

% create folder to store data in
% -> to data\sys#\ref0_<>\dRe\<subdir1>\<subdir2>\seed_<seed>.mat
str_iRe = iRe2str(iRe,nRe);
subdir2 = sprintf('%s_Re_%.2e',str_iRe,Re);
subdir2 = replace(subdir2,'.','p'); % replace . with p
subdir2 = fullfile(subdir1,subdir2);
if ~isfolder(subdir2)
    mkdir(subdir2);
end

%% run simulations
[s.opts, s.u0, s.y0, s.xcl0, s.Z, s.Lf, s.Cz, s.Tcl, s.u_cl, s.y_cl] ...
    = main_MC(opts,sigs,plant,Cz0,Tcl0,yr0,s.e0,yr1,ur1,s.e1);

%% save data
% ----------------------- save data for run iRe & seed ----------------
fn = sprintf('seed_%d.mat',seed);
fn = fullfile(subdir2,fn);
fn_short = strrep(fn,proj_dir,'');
fprintf('Saving data to file: \n\t%s \n',fn_short);
% save(fn,'opts','e0','e1','u0','y0','xcl0',...
%     'Z','Lf','Cz','Tcl','u_cl','y_cl');
save(fn,"-fromstruct",s)
fprintf('File saved successfully!\n');
end

% get zero-padded iRe
function str_iRe = iRe2str(iRe,nRe)
nDigits = ceil(log10(nRe + 1));
format = ['%0', num2str(nDigits), 'd'];
str_iRe = sprintf(format,iRe);
end