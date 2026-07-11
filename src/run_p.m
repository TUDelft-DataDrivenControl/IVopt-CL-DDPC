function run_p(ii,opts2,spP,nP,seeds,p_all,f,N,Ncl,ny,nu,Re,plant,subdir1,sigs,Cz0,Tcl0,proj_dir,yr0_full)
%% This function acts as a wrapper for main_MC.m for a single value of p

sP = struct;
sP.opts = opts2;

% ------------------------------ define p run -----------------------------
[~,iP] = ind2sub([spP,nP],ii);
p = p_all(iP);
sP.opts.p = p;

% ----- initial CL-sim length & reference -----------
[Nbar, sP.Nbar] = deal(p + f + (N-1)*opts2.DMCS); % sim. length of initial controller
% select initial reference (yr0) from pre-generated full reference
sP.yr0 = yr0_full(:,1:Nbar);

% ---------- references for subsequent closed-loop simulations ------------
sP.yr1 = make_reference(Ncl+f,ny); % y-ref
P0  = dcgain(plant(:,1:nu));       % DC gain
sP.ur1 = P0\sP.yr1;                % u-ref

% ----------------------- save settings for run iP ------------------------
% -> to data\sys#\ref#\dp\<subdir1>\<subdir2>\<iP>_settings.mat
str_iP = iN2str(iP,nP); % zero-padded <iP> based on # of decimals for nP
subdir2 = sprintf('%s_p_%d',str_iP,p); % subdir2 name
subdir2 = fullfile(subdir1,subdir2);  % subdir2 path
if ~isfolder(subdir2)
    mkdir(subdir2); % create subdir2
end
nset_sd2 = append(str_iP,'_settings.mat');
if ~isfile(fullfile(subdir2,nset_sd2))
    save(fullfile(subdir2,nset_sd2),"-fromstruct",sP) % 'opts','Nbar','yr0','yr1','ur1'
end

% ------------------------- define seed run & noise -----------------------
sE = struct('opts',sP.opts);
seed = seeds(ii); sE.opts.seed = seed;
rng(seed);
sE.e0 = mvnrnd(zeros(ny,1),Re,Nbar).'; % innovation noise
sE.e1 = mvnrnd(zeros(ny,1),Re,Ncl).';  % innovation noise

%% run simulations
[sE.opts, sE.u0, sE.y0, sE.xcl0, sE.Z, sE.Lf, sE.Cz, sE.Tcl, sE.u_cl, sE.y_cl, sE.Cases] ...
    = main_MC(sE.opts,sigs,plant,Cz0,Tcl0,sP.yr0,sE.e0,sP.yr1,sP.ur1,sE.e1);

%% save data
fn = sprintf('seed_%d.mat',seed);
fn = fullfile(subdir2,fn);
fn_short = strrep(fn,proj_dir,'');
fprintf('Saving data to file: \n\t%s \n',fn_short);
save(fn,"-fromstruct",sE);
%'opts','e0','e1','u0','y0','xcl0',...
% 'Z','Lf','Cz','Tcl','u_cl','y_cl','Cases'
fprintf('File saved successfully!\n');
end

% get zero-padded iN
function str_iN = iN2str(iN,nN)
nDigits = ceil(log10(nN + 1));
format = ['%0', num2str(nDigits), 'd'];
str_iN = sprintf(format,iN);
end