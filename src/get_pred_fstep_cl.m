function [sys_cl,x0_cl] = get_pred_fstep_cl(f,p,ny,nu,ICuf,ICyf,ICup,ICyp,ICyr,ICur,Tf_u,Gf_tKp_u,Gf_tKp_y,up,yp,yrf,urf)
% get closed-loop sysem with yf-predictor and controller

Uf2Uf0 = ICuf + ICyf*Tf_u;
left_Uf = eye(nu*f) - Uf2Uf0;
inv_Uf = left_Uf\[ICup+ICyf*Gf_tKp_u ICyp+ICyf*Gf_tKp_y ICyr ICur];
C2up  = inv_Uf(:,1:size(ICup,2));
C2yp  = inv_Uf(:,size(ICup,2)+1:size(ICup,2)+size(ICyp,2));
MIyr  = inv_Uf(:,end-size([ICyr ICur],2)+1:end-size(ICur,2));
MIur  = inv_Uf(:,end-size(ICur,2)+1:end);

Tup = Gf_tKp_u + Tf_u*C2up;
Typ = Gf_tKp_y + Tf_u*C2yp;
Tyr = Tf_u*MIyr;
Tur = Tf_u*MIur;

% creating the A matrix
nxL = (p+2*f-1)*(nu+ny);
Aup = [zeros(p*nu,f*nu) eye(p*nu)]*[eye(p*nu,nxL);C2up C2yp MIur MIyr];
Ayp = [zeros(p*ny,f*ny) eye(p*ny)]*[zeros(p*ny,p*nu) eye(p*ny) zeros(p*ny,(nu+ny)*(2*f-1)); Tup Typ Tur Tyr];
ABur = [zeros((2*f-1)*nu,f*nu) eye((2*f-1)*nu)];
Aur = ABur(:,1:(2*f-1)*nu); Aur = [zeros((2*f-1)*nu,p*(nu+ny)) Aur zeros((2*f-1)*nu,(2*f-1)*ny)];
AByr = [zeros((2*f-1)*ny,f*ny) eye((2*f-1)*ny)];
Ayr = AByr(:,1:(2*f-1)*ny); Ayr = [zeros((2*f-1)*ny,nxL-(2*f-1)*ny) Ayr];
%     Aur = [zeros((2*f-2)*nu,p*(nu+ny)+nu) eye((2*f-2)*nu) zeros((2*f-2)*nu,(2*f-1)*ny);zeros(nu,nxL)];
%     Ayr = [zeros((2*f-2)*ny,p*(nu+ny)+(2*f-1)*nu+ny) eye((2*f-2)*ny);zeros(ny,nxL)];
A = [Aup;Ayp;Aur;Ayr];

% creating the B matrix
B = zeros(nxL,(nu+ny)*f);
B(p*(nu+ny)+1:p*(nu+ny)+(2*f-1)*nu,1:nu*f) = ABur(:,end-nu*f+1:end);
B(p*(nu+ny)+(2*f-1)*nu+1:end,nu*f+1:end) = AByr(:,end-nu*f+1:end);

% creating the C matrix
C = [C2up C2yp MIur MIyr; Tup Typ Tur Tyr];

% creating the D matrix
D = zeros((nu+ny)*f);

% creating closed-loop system
sys_cl = ss(A,B,C,D,-1);
x0_cl = [up(:);yp(:);urf(:);yrf(:)]; % urf: nu x 2*f-1, yrf: ny x 2*f-1

%% ensure stability of closed-loop system
eigs = pole(sys_cl);
itk = 0;
while any(abs(eigs) > 1)
    itk = itk + 1;
    if itk == 1
        % transform system to eigenvalue-revealing modal form
        [sys_cl2,T] = canon(sys_cl,'modal');

        % identify blocks on diagonal
        Dblks = identify_diag_blocks(sys_cl2.A);
        
        % pre-define for for-loop over blocks on diagonal:
        eigs = cell(size(Dblks));       % eigenvalues
        abs_eigs = nan(size(Dblks));    % absolute value of eigenvalues
        angle_eigs = eigs;              % angle of eigenvalues in complex plane
        sblkJ = nan(size(Dblks));       % size of blocks on diagonal
        for kk = 1:length(Dblks) % loop over blocks on diagonal
            blkJ = Dblks{kk};               % define block
            sblkJ(kk) = size(blkJ,1);       % get block size
            
            if sblkJ(kk) == 1
                % if block is of size 1 -> it's an eigenvalue
                eigs{kk}=blkJ;
                abs_eigs(kk) = abs(blkJ);
            else
                eigs{kk} = unique(eig(blkJ)); % get unique eigenvalues
                unique_abs_eigs = unique(abs(round(eigs{kk}, 8))); % get unique abs. val. of eigenvalues

                % there should be only one unique absolute eigenvalue, but
                % if it's less than 1, it's not going to cause any issues
                throw_error = max(unique_abs_eigs,[],'all') > 1 && numel(unique_abs_eigs) > 1;
                if throw_error
                    error('block encountered with non-distinct eigenvalues of which at least one with an absolute eigenvalue greater than one')
                elseif numel(unique_abs_eigs) > 1
                    % select the maximum eigenvalue
                    abs_eigs(kk) = max(unique_abs_eigs,[],'all');
                else
                    abs_eigs(kk) = unique_abs_eigs;
                end
            end
            angle_eigs{kk} = angle(eigs{kk});
        end
    end
    msk = find(abs_eigs > 1); % find eigenvalues that must be reflected
    abs_eigs_old = abs_eigs;
    abs_eigs(msk) = 1-rem(abs_eigs(msk),1);
    abs_eigs(msk);
    for kk = 1:length(msk)
        msk_idx = msk(kk);
        % adjust blocks on diagonal
        Dblks{msk_idx} = Dblks{msk_idx}*abs_eigs(msk_idx)/abs_eigs_old(msk_idx);
    end
    Anew = T\blkdiag(Dblks{:})*T;
    eigs = eig(Anew);
end

% adjust C matrix if A has been adjusted
if itk > 0
    % determine how many last block-rows are affected
    nrows = min(f,p);
    C(nu*f-nrows*nu+1:nu*f,:) = Anew(nu*p-nrows*nu+1:nu*p,:);
    C(end-nrows*ny+1:end,:)   = Anew((nu+ny)*p-nrows*ny+1:(nu+ny)*p,:);
    sys_cl = ss(Anew,B,C,D,-1);
end

end

%% Helper functions
function Dblocks = identify_diag_blocks(A)
    % Identify Jordan blocks in A, even if they are not strictly adjacent.
    
    Dblocks = {};  % Cell array to store Jordan blocks
    
    % Iterate over all rows to find the diagonal block(s)
    B = 1; % block begin
    E = B; % block end
    R = size(A,1);
    while B<=R
        while E<R && any(any(A(E+1:R,B:E)|A(B:E,E+1:R).'))
            E = E+1;
        end
        Dblocks{end+1} = A(B:E,B:E); %%#ok<SAGROW>
        E = E+1;
        B = E;
    end
    assert(all(size(blkdiag(Dblocks{:}))==size(A)),'Jblocks do not correspond to A matrix');
%     % check if eigenvalues are unique
%     allunique = false;
%     while ~allunique
%         for k = 1:length(Jblocks)
%             Jblock = Jblocks{k};
%             uniqueAbsEigVals = unique(abs(round(eig(Jblock), 8)));  % Avoid numerical errors
%             if any(uniqueAbsEigVals > 1) % otherwise it won't create problems
%                 if (istriu(Jblock) || istril(Jblock)) && numel(uniqueAbsEigVals) > 1
%                     Jblocks2 = Jblocks; Jblocks2{end+1} = {};
%                     Jblocks2(1:k-1) = Jblocks(1:k-1);
%                     Jblocks2(k+2:end) = Jblocks(k+1:end);
%                     dJblock = diag(Jblock);
%                     i_new = find(dJblock~=dJblock(1),1,'first');
%                     Jblocks2{k} = Jblock(1:i_new-1,1:i_new-1);
%                     Jblocks2{k+1} = Jblock(i_new:end,i_new:end);
%                     Jblocks = Jblocks2; clear Jblocks2;
%                     break;
%                 end
%             end
%             if k==length(Jblocks)
%                 allunique = true;
%             end
%         end
%     end

end