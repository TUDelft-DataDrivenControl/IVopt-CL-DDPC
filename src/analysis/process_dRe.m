%% iterate over Re values - initial data processing
subdir2s = dir(pwd);
isub = [subdir2s(:).isdir]; 
subdir2s = {subdir2s(isub).name};
subdir2s = subdir2s(~ismember(subdir2s,{'.','..','mfiles'}));

nRe = numel(subdir2s);
nX = nRe;

%% initializing measures
Uf_ivs = {'iv1','iv2a','iv2c','iv3c','iv4a','iv4c','iv5a','iv5c','iv6c'};
Yf_ivs = {'iv2b','iv3a','iv4b','iv5b','iv6a'};
num_Uf_ivs = numel(Uf_ivs);
num_Yf_ivs = numel(Yf_ivs);

m0_Uf_mean = cell(num_Uf_ivs,nX);
m0_Yf_mean = cell(num_Yf_ivs,nX);
m0_Uf_median = m0_Uf_mean;
m0_Yf_median = m0_Yf_mean;
m0_Uf_pctiles = cell(num_Uf_ivs,nX);
m0_Yf_pctiles = cell(num_Yf_ivs,nX);

pctiles = 0:5:100;
num_pctiles  = numel(pctiles);
m1_Uf_data     = zeros(num_Uf_ivs, nX, spX);
m1_Yf_data     = zeros(num_Yf_ivs, nX, spX);
iyf = nu*f + (1:ny*f);

Cases = {'iv1','iv2a','iv2b','iv2c','iv3a','iv3c','iv4a','iv4b','iv4c', ...
         'iv5a','iv5b','iv5c','iv6a','iv6c','CLSPC','actLf'};
num_Cases = numel(Cases);
IDerrorTypes = {'Up','Yp','Uf'};
m2_data    = zeros(numel(IDerrorTypes), numel(Cases), nX, spX);
num_IDerrorTypes = numel(IDerrorTypes);
cost_types = {'cost_u','cost_y','cost_tot'};
nCostTypes = numel(cost_types);
m3_data    = zeros(nCostTypes, numel(Cases), nX, spX);

%% loop over Re values
for iRe = 1:nRe
    Re = Re_all(iRe);
    iX = iRe;
    subdir2 = choose_subdir_by_iX(subdir2s, iX);
    cd(subdir2);
    num_diags = f+N-1;
    clc;
    fprintf('Processing Re index %d/%d (Re = %g) in subdir: %s\n', iRe, nRe, Re, subdir2);
    parfor kIVu = 1:num_Uf_ivs % parfor
        iv_name = Uf_ivs{kIVu};
        m0_Uf_mean2 = zeros(nu,num_diags);
        m0_Uf_median2 = m0_Uf_mean2;
        m0_Uf_pctiles2 = zeros(nu,num_diags,num_pctiles);
        m0_Uf_data2 = zeros(nu*f,N,spX);
        for ks = 1:spRe
            seed = seeds(ks,iRe);
            fndata = sprintf('seed_%d.mat',seed);
            Z = load(fndata,'Z').Z;
            Uf_iv = Z.([iv_name,'_']);
            m0_Uf_data2(:,:,ks) = Uf_iv;
        end
        ij_adiags = get_subind_diags(f,N,nr=nu,anti=true);
        for kd = 1:num_diags
            rows = ij_adiags{kd}(:,1);
            cols = ij_adiags{kd}(:,2);
            rows2 = repmat(rows,spX,1);
            cols2 = repmat(cols,spX,1);
            i3s = kron( (1:spX).', ones(numel(cols),1) );
            idxlin = sub2ind(size(m0_Uf_data2),rows2,cols2,i3s);
            Uf_sel = m0_Uf_data2(idxlin);
            Uf_sel = reshape(Uf_sel,nu,[]);
            m0_Uf_mean2(:,kd)      = mean(Uf_sel,2);
            m0_Uf_median2(:,kd)    = median(Uf_sel,2);
            m0_Uf_pctiles2(:,kd,:) = prctile(Uf_sel,pctiles,2);
        end
        m0_Uf_mean{kIVu,iX}    = m0_Uf_mean2;
        m0_Uf_median{kIVu,iX}  = m0_Uf_median2;
        m0_Uf_pctiles{kIVu,iX} = m0_Uf_pctiles2;
        
        fprintf('  [Re %d/%d] Uf IV %d/%d done\n', iRe, nRe, kIVu, num_Uf_ivs);
    end
    fprintf('  [Re %d/%d] All Uf IVs done\n', iRe, nRe);
    parfor kIVy = 1:num_Yf_ivs % parfor
        iv_name = Yf_ivs{kIVy};
        m0_Yf_mean2 = zeros(ny,num_diags);
        m0_Yf_median2 = m0_Yf_mean2;
        m0_Yf_pctiles2 = zeros(ny,num_diags,num_pctiles);
        m0_Yf_data2 = zeros(ny*f,N,spX);
        for ks = 1:spRe
            seed = seeds(ks,iRe);
            fndata = sprintf('seed_%d.mat',seed);
            Z = load(fndata,'Z').Z;
            switch iv_name
                case 'iv3a'
                    Yf_iv = Z.iv3a_(1:ny*f,:);
                case 'iv6a'
                    Yf_iv = Z.iv6a_;
                otherwise
                    Yf_iv = Z.([iv_name,'_'])(iyf,:);
            end
            m0_Yf_data2(:,:,ks) = Yf_iv;
        end
        ij_adiags = get_subind_diags(f,N,nr=ny,anti=true);
        for kd = 1:num_diags
            rows = ij_adiags{kd}(:,1);
            cols = ij_adiags{kd}(:,2);
            rows2 = repmat(rows,spX,1);
            cols2 = repmat(cols,spX,1);
            i3s = kron( (1:spX).', ones(numel(cols),1) );
            idxlin = sub2ind(size(m0_Yf_data2),rows2,cols2,i3s);
            Yf_sel = m0_Yf_data2(idxlin);
            Yf_sel = reshape(Yf_sel,ny,[]);
            m0_Yf_mean2(:,kd)      = mean(Yf_sel,2);
            m0_Yf_median2(:,kd)    = median(Yf_sel,2);
            m0_Yf_pctiles2(:,kd,:) = prctile(Yf_sel,pctiles,2);
        end
        m0_Yf_mean{kIVy,iX}    = m0_Yf_mean2;
        m0_Yf_median{kIVy,iX}  = m0_Yf_median2;
        m0_Yf_pctiles{kIVy,iX} = m0_Yf_pctiles2;

        fprintf('  [Re %d/%d] Yf IV %d/%d done\n', iRe, nRe, kIVy, num_Yf_ivs);
    end
    fprintf('  [Re %d/%d] All Yf IVs done\n', iRe, nRe);
    parfor ks = 1:spRe % parfor
        seed = seeds(ks,iRe);
        fndata = sprintf('seed_%d.mat',seed);
        [Cases,Cz,FroIDerror,Lf,Tcl,Z,cost_tot,cost_u,cost_u1,cost_u2,cost_y,...
          e0,e1,opts,u0,u_cl,u_iv,xcl0,y0,y_cl,y_iv] = load_seedmat(fndata);
        for kIVu = 1:num_Uf_ivs
            iv_name = Uf_ivs{kIVu};
            switch iv_name
                case 'iv2a'
                otherwise
                    Uf_iv = Z.([iv_name,'_']);
                    m1_Uf_data(kIVu, iX, ks) = norm(Uf_iv - Z.iv2a_, 'fro');
            end
        end
        for kIVy = 1:num_Yf_ivs
            iv_name = Yf_ivs{kIVy};
            switch iv_name
                case 'iv2b'
                    Yf_iv = 0;
                    calc_norm = false;
                case 'iv3a'
                    Yf_iv = Z.iv3a_(1:ny*f,:);
                    calc_norm = true;
                case 'iv6a'
                    Yf_iv = Z.iv6a_;
                    calc_norm = true;
                otherwise
                    Yf_iv = Z.([iv_name,'_'])(iyf,:);
                    calc_norm = true;
            end
            if calc_norm
                m1_Yf_data(kIVy, iX, ks) = norm(Yf_iv - Z.iv2b_(iyf,:), 'fro');
            end
        end
        for kType = 1:num_IDerrorTypes
            IDerrorType = IDerrorTypes{kType};
            for kIVn = 1:num_Cases
                IVn = Cases{kIVn};
                m2_data(kType, kIVn, iX, ks) = FroIDerror.(IVn).(IDerrorType);
            end
        end
        for kType = 1:nCostTypes
            cost_type = cost_types{kType};
            switch cost_type
                case 'cost_u'
                    cost = cost_u;
                case 'cost_y'
                    cost = cost_y;
                otherwise
                    cost = cost_tot;
            end
            for kIVn = 1:num_Cases
                IVn = Cases{kIVn};
                m3_data(kType, kIVn, iX, ks) = cost.(IVn);
            end
        end
        if mod(ks, max(1, floor(spRe/10))) == 0 || ks == spRe
            fprintf('  [Re %d/%d] Seed %d/%d done\n', iRe, nRe, ks, spRe);
        end
    end
    fprintf('  [Re %d/%d] All seeds done\n', iRe, nRe);
    cd(subdir1);
end

fprintf("Processing m0 data\n")
m0 = struct;
for k = 1:numel(Uf_ivs)
    for iX = 1:nX
        iXstr = sprintf('iX%d',iX);
        m0.Uf.(Uf_ivs{k}).(iXstr).mean     = m0_Uf_mean{k,iX};
        m0.Uf.(Uf_ivs{k}).(iXstr).median   = m0_Uf_median{k,iX};
        m0.Uf.(Uf_ivs{k}).(iXstr).pctiles  = m0_Uf_pctiles{k,iX};
    end
end
for k = 1:numel(Yf_ivs)
    for iX = 1:nX
        iXstr = sprintf('iX%d',iX);
        m0.Yf.(Yf_ivs{k}).(iXstr).mean     = m0_Yf_mean{k,iX};
        m0.Yf.(Yf_ivs{k}).(iXstr).median   = m0_Yf_median{k,iX};
        m0.Yf.(Yf_ivs{k}).(iXstr).pctiles  = m0_Yf_pctiles{k,iX};
    end
end

fprintf("Processing m1 data\n")
m1 = struct();
for k = 1:numel(Uf_ivs)
    m1.Uf.(Uf_ivs{k}).data     = squeeze( m1_Uf_data(k,:,:) );
    m1.Uf.(Uf_ivs{k}).mean     = mean(    m1.Uf.(Uf_ivs{k}).data, 2);
    m1.Uf.(Uf_ivs{k}).median   = median(  m1.Uf.(Uf_ivs{k}).data, 2);
    m1.Uf.(Uf_ivs{k}).pctiles  = prctile( m1.Uf.(Uf_ivs{k}).data, pctiles, 2);
end
for k = 1:numel(Yf_ivs)
    m1.Yf.(Yf_ivs{k}).data     = squeeze( m1_Yf_data(k,:,:) );
    m1.Yf.(Yf_ivs{k}).mean     = mean(    m1.Yf.(Yf_ivs{k}).data, 2);
    m1.Yf.(Yf_ivs{k}).median   = median(  m1.Yf.(Yf_ivs{k}).data, 2);
    m1.Yf.(Yf_ivs{k}).pctiles  = prctile( m1.Yf.(Yf_ivs{k}).data, pctiles, 2);
end

fprintf("Processing m2 data\n")
m2 = struct();
for kType = 1:numel(IDerrorTypes)
    typeName = IDerrorTypes{kType};
    for kIVn = 1:numel(Cases)
        caseName = Cases{kIVn};
        m2.(typeName).(caseName).data    = squeeze( m2_data(kType, kIVn, :, :) );
        m2.(typeName).(caseName).mean    = mean(    m2.(typeName).(caseName).data, 2);
        m2.(typeName).(caseName).median  = median(  m2.(typeName).(caseName).data, 2);
        m2.(typeName).(caseName).pctiles = prctile( m2.(typeName).(caseName).data, pctiles, 2);
    end
end

fprintf("Processing m3 data\n")
m3 = struct();
for kType = 1:nCostTypes
    costName = cost_types{kType};
    for kIVn = 1:numel(Cases)
        caseName = Cases{kIVn};
        m3.(costName).(caseName).data    = squeeze( m3_data(kType, kIVn, :, :) );
        m3.(costName).(caseName).mean    = mean(    m3.(costName).(caseName).data, 2);
        m3.(costName).(caseName).median  = median(  m3.(costName).(caseName).data, 2);
        m3.(costName).(caseName).pctiles = prctile( m3.(costName).(caseName).data, pctiles, 2);
    end
end

fndata = 'processed_data.mat';
fprintf('Saving data to %s\n',fndata);
save(fndata,'m0','m1','m2','m3');


function chosenDir = choose_subdir_by_iX(subdirs, iX)
% CHOOSE_SUBDIR_BY_IX selects the subdirectory whose leading zero-padded number matches iX.
% Example: for iX=3, matches '003_someName' if present.
    chosenDir = '';
    for i = 1:numel(subdirs)
        dirName = subdirs{i};
        tokens = regexp(dirName, '^(\d+)', 'tokens');
        if ~isempty(tokens)
            dirNum = str2double(tokens{1}{1});
            if dirNum == iX
                chosenDir = dirName;
                return;
            end
        end
    end
    if isempty(chosenDir)
        warning('No subdirectory found starting with number %d.', iX);
    end
end

function [Cases,Cz,FroIDerror,Lf,Tcl,Z,cost_tot,cost_u,cost_u1,cost_u2,cost_y,...
          e0,e1,opts,u0,u_cl,u_iv,xcl0,y0,y_cl,y_iv] = load_seedmat(fnpath)
    load(fnpath);
end
