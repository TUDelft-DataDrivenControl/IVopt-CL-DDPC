function make_solver(obj,usr_con)
% makes solver based on the provided constraints, provided as the following
% fields of the usr_con structure:
% - u_min:  minimum inputs          nu x 1 (double)
% - u_max:  maximum inputs          nu x 1 (double)
% - y_min:  minimum outputs         ny x 1 (double)
% - y_max:  maximum outputs         ny x 1 (double)
% - du_max: maximum |u_k-u_{k-1}|   nu x 1 (double)
% - dy_max: maximum |y_k-y_{k-1}|   ny x 1 (double)
% - expr: cell array of casadi.SX expressions that indicate other
%         user-defined constraints. The employed variables (uf, yf, rf, u0, y0)
%         must be provided if used (see below). Either in form F1. or F2.
%         F1:   nexpr x 2 (cell array), where each row is of the form:
%               {expr, relationship w.r.t. 0 ('<=','==',or '>=')}
%         F2: 2*nexpr x 1 (cell array). vectorized form 1
% - uf:     future inputs           nu x f (casadi.SX)
% - yf:     future outputs          ny x f (casadi.SX)
% - rf:     future reference        ny x f (casadi.SX)
% - u0:     last past input         nu x 1 (casadi.SX)
% - y0:     last past output        ny x 1 (casadi.SX)

% specify solver method
obj.solve = @obj.optimizer_solve;

%% Preliminaries
% possible constraint structure field types
str1 = {'uf','yf','rf','u0','y0'}; %-> if present, so should be expr field
str2 = {'u_min','u_max','y_min','y_max','du_max','dy_max'};

% initializing upper, lower bounds, constraint matrix
lbx = []; ubx = []; A = []; lba = []; uba = [];

%% parsing user-defined constraints
usr_con = parse_usr_con(obj, usr_con,str1,str2);

%% make variables/parameters up, yp, rf, uf, yf

if isfield(usr_con,'u0')
    obj.Prob.up_ = [casadi.SX.sym('up_endmin1', obj.nu, obj.p-1),usr_con.u0];
else
    obj.Prob.up_ = casadi.SX.sym('up', obj.nu, obj.p);
end
if isfield(usr_con,'y0')
    obj.Prob.yp_ = [casadi.SX.sym('yp_endmin1', obj.ny, obj.p-1),usr_con.y0];
else
    obj.Prob.yp_ = casadi.SX.sym('yp', obj.ny, obj.p);
end
if isfield(usr_con,'rf')
    obj.Prob.rf_ = usr_con.rf;
else
    obj.Prob.rf_ = casadi.SX.sym('rf', obj.ny, obj.f);
end
if isfield(usr_con,'uf')
    obj.Prob.uf_ = usr_con.uf;
else
    obj.Prob.uf_ = casadi.SX.sym('uf', obj.nu, obj.f);
end
if isfield(usr_con,'yf')
    obj.Prob.yf_ = usr_con.yf;
else
    obj.Prob.yf_ = casadi.SX.sym('yf', obj.ny, obj.f);
end

%% parameters - up, yp, rf, Lu, Ly, Gu -> p
if obj.options.ExplicitPredictor
    obj.Prob.Lu_ = casadi.SX.sym('Lu', obj.f*obj.ny, obj.p*obj.nu);
    obj.Prob.Ly_ = casadi.SX.sym('Ly', obj.f*obj.ny, obj.p*obj.ny);
    obj.Prob.Gu_ = casadi.SX.sym('Gu', obj.f*obj.ny, obj.f*obj.nu);
    obj.Prob.p_ = [obj.Prob.up_(:); obj.Prob.yp_(:); obj.Prob.rf_(:); obj.Prob.Lu_(:); obj.Prob.Ly_(:); obj.Prob.Gu_(:)];
    obj.Prob.get_p = casadi.Function('get_p',... function to optain parameters
        {obj.Prob.up_,obj.Prob.yp_,obj.Prob.rf_,obj.Prob.Lu_,obj.Prob.Ly_,obj.Prob.Gu_},... parameters
        {obj.Prob.p_},... vector with parameters
        {'up','yp','rf','Lu','Ly','Gu'},{'p'}); % naming
    obj.Prob.p2Gu = casadi.Function('p2Gu',{obj.Prob.p_},{obj.Prob.Gu_});
    obj.Prob.p2Ly = casadi.Function('p2Ly',{obj.Prob.p_},{obj.Prob.Ly_});
    obj.Prob.p2Lu = casadi.Function('p2Lu',{obj.Prob.p_},{obj.Prob.Lu_});
else
    if obj.options.use_IV
        m1 = obj.pfid*obj.nu + obj.p*obj.ny;
        m2 = m1 + obj.fid*obj.ny;
    else
        m1 = obj.N;
        m2 = obj.pfid*(obj.nu+obj.ny);
    end
    
    obj.Prob.LHS_ = casadi.SX.sym('LHS', m2, m1);
    obj.Prob.p_ = [obj.Prob.up_(:); obj.Prob.yp_(:); obj.Prob.rf_(:); obj.Prob.LHS_(:)];
    obj.Prob.get_p = casadi.Function('get_p',... function to optain parameters
        {obj.Prob.up_,obj.Prob.yp_,obj.Prob.rf_,obj.Prob.LHS_},... parameters
        {obj.Prob.p_},... vector with parameters
        {'up','yp','rf','LHS'},{'p'}); % naming
end
zero_p = zeros(size(obj.Prob.p_));

%% optimization variables - uf, yf, G -> x
if obj.options.ExplicitPredictor
    obj.Prob.x_ = [obj.Prob.uf_(:); obj.Prob.yf_(:)];
else
    obj.Prob.G_ = casadi.SX.sym('G', m1, obj.nGcols);
    obj.Prob.x_ = [obj.Prob.uf_(:); obj.Prob.yf_(:); obj.Prob.G_(:)];
    
    % make functions to get results
    obj.Prob.x2G = @(x) reshape(x((obj.nu+obj.ny)*obj.f+1:(obj.nu+obj.ny)*obj.f+m1*obj.nGcols),m1,obj.nGcols);
end
obj.Prob.x2uf = @(x) reshape(x(1:obj.nu*obj.f),obj.nu,obj.f); % define in terms of beginning s.t. useable for result from QP3
obj.Prob.x2yf = @(x) reshape(x(obj.nu*obj.f+1:(obj.nu+obj.ny)*obj.f),obj.ny,obj.f);
zero_x = zeros(size(obj.Prob.x_));

%% process usr_con.expr
if isfield(usr_con,'expr')
    % obtain con_LHS & con_gleq from usr_con.expr
    if isvector(usr_con.expr)
        con_LHS  = usr_con.expr(1:2:end); con_LHS = con_LHS(:);
        con_gleq = usr_con.expr(2:2:end); con_gleq= con_gleq(:);
    elseif ismatrix(usr_con.expr)
        con_LHS  = usr_con.expr(:,1);
        con_gleq = usr_con.expr(:,2);
    end
    
    % initialize lbx, ubx
    for k = 1:numel(con_LHS)
        lhs = con_LHS{k}(:);

        Ax_new = sparse(casadi.DM(jacobian(lhs,obj.Prob.x_)));
        % simple check for identity matrix: constraint of form lbx <= x <= ubx
        [i,j,a] = find(Ax_new);
        i_counts = grouptransform(i,i,@numel);
        bx_msk = i_counts == 1 & a == 1;
        i_bx = i(bx_msk);
        j_bx = j(bx_msk);
        i_ax = unique(i(~bx_msk));
        
        Ap_new = casadi.DM(jacobian(lhs,obj.Prob.p_));
        const_new = casadi.substitute(lhs,obj.Prob.x_,zero_x);
        const_new = casadi.substitute(const_new,obj.Prob.p_,zero_p);
        p_terms   = Ap_new*obj.Prob.p_;
        rhs       = -p_terms-const_new;
        num_rows  = size(const_new,1);
        
        % set new upper/lower bounds
        if strcmp(con_gleq{k},'==')
            lb_new = rhs;
            ub_new = rhs;
        elseif strcmp(con_gleq{k},'<=') || strcmp(con_gleq{k},'=<')
            lb_new = -inf(num_rows,1);
            ub_new = rhs;
        elseif strcmp(con_gleq{k},'>=') || strcmp(con_gleq{k},'=>')
            lb_new = rhs;
            ub_new = inf(num_rows,1);
        else
            error('Constraint specified incorrectly.')
        end
        
        % constraints of the form  lba <= Ax <= uba
        if ~isempty(i_ax)
            A   = [A;  Ax_new(i_ax,:)];
            lba = [lba;lb_new(i_ax,1)];
            uba = [uba;ub_new(i_ax,1)];
        end
        
        % constraints of the form lbx <= x <= ubx
        if ~isempty(i_bx)
            ubx_new = casadi.SX.inf(size(obj.Prob.x_)); ubx_new(j_bx,1) = ub_new(i_bx,1);
            lbx_new = -ubx_new;                         lbx_new(j_bx,1) = lb_new(i_bx,1);
            lbx = [lbx lbx_new]; % will take max(<>,[],2)
            ubx = [ubx ubx_new]; % will take min(<>,[],2)
        end
    end
end
%% process usr_con str2 fields
% add missing fields
str2 = {'u_min','u_max','y_min','y_max','du_max','dy_max'};
for k=1:length(str2)
    fn = str2{k};
    if ~isfield(usr_con,fn)
        switch fn
            case 'u_min'
                usr_con.u_min = -inf(obj.nu,1);
            case 'u_max'
                usr_con.u_max =  inf(obj.nu,1);
            case 'y_min'
                usr_con.y_min = -inf(obj.ny,1);
            case 'y_max'
                usr_con.y_max =  inf(obj.ny,1);
            case {'du_max','dy_max'}
                usr_con.(fn) = [];
        end
    end
end
if obj.options.ExplicitPredictor
    lbx = [lbx, [repmat(usr_con.u_min,obj.f,1); repmat(usr_con.y_min,obj.f,1)]];
    ubx = [ubx, [repmat(usr_con.u_max,obj.f,1); repmat(usr_con.y_max,obj.f,1)]];
else
    lbx = [lbx, [repmat(usr_con.u_min,obj.f,1); repmat(usr_con.y_min,obj.f,1);-inf(numel(obj.Prob.G_),1)]];
    ubx = [ubx, [repmat(usr_con.u_max,obj.f,1); repmat(usr_con.y_max,obj.f,1); inf(numel(obj.Prob.G_),1)]];
end
if ~isempty(usr_con.du_max)
    lba = [lba; -usr_con.du_max+obj.Prob.up_(:,end); -repmat(usr_con.du_max,obj.f-1,1)];
    uba = [uba;  usr_con.du_max+obj.Prob.up_(:,end);  repmat(usr_con.du_max,obj.f-1,1)];
    du_ = [obj.Prob.uf_(:,1)-obj.Prob.up_(:,end) obj.Prob.uf_(:,2:end)-obj.Prob.uf_(:,1:end-1)];
    A   = [A;casadi.DM(jacobian(du_(:),obj.Prob.x_))];
end
if ~isempty(usr_con.dy_max)
    lba = [lba;-usr_con.dy_max+obj.Prob.yp_(:,end); -repmat(usr_con.dy_max,obj.f-1,1)];
    uba = [uba; usr_con.dy_max+obj.Prob.yp_(:,end);  repmat(usr_con.dy_max,obj.f-1,1)];
    dy_ = [obj.Prob.yf_(:,1)-obj.Prob.yp_(:,end) obj.Prob.yf_(:,2:end)-obj.Prob.yf_(:,1:end-1)];
    A   = [A;casadi.DM(jacobian(dy_(:),obj.Prob.x_))];
end

%% remove yf from optimization variable -> integrate into cost, constraints
er_ = obj.Prob.yf_ - obj.Prob.rf_; % error w.r.t. reference
du_ = horzcat(obj.Prob.uf_(:,1)    -obj.Prob.up_(:,end), ...
              obj.Prob.uf_(:,2:end)-obj.Prob.uf_(:,1:end-1)); % u_{k+1}-u_k
obj.Prob.cost =   er_(:).'*obj.Prob.Q *er_(:) ...
       + obj.Prob.uf_(:).'*obj.Prob.R *obj.Prob.uf_(:) ...
                + du_(:).'*obj.Prob.dR*du_(:);

%% Constraints - Dynamics
if obj.options.ExplicitPredictor
    ulba_dyn = obj.Prob.Lu_*obj.Prob.up_(:)+obj.Prob.Ly_*obj.Prob.yp_(:);
    Adyn = [-obj.Prob.Gu_ speye(obj.ny*obj.f)];
else
    Hf_= [obj.make_CasADi_Hankel([obj.Prob.up_ obj.Prob.uf_],obj.pfid,obj.nGcols,'u');...
          obj.make_CasADi_Hankel([obj.Prob.yp_ obj.Prob.yf_],obj.pfid,obj.nGcols,'y')];
    x_small = [obj.Prob.uf_(:);obj.Prob.yf_(:)]; % x_ without G
    Adyn = [casadi.DM(-jacobian(Hf_(:),x_small)), kron(speye(obj.nGcols),obj.Prob.LHS_)];
    ulba_dyn = casadi.substitute(Hf_(:),x_small,zeros(size(x_small)));
    
    obj.Prob.get_Adyn = casadi.Function('get_Adyn',{obj.Prob.p_},{Adyn}); % needed for initial guess in optimizer_solve
end

%% make functions for regular solver
% solver options
opts = obj.options.opts;
if isempty(opts)
    % leave obj.Prob.cas_opts as default
elseif isstruct(opts)
    if ~isfield(opts,'solver')
        opts.solver = obj.Prob.cas_opts.solver;
    end
    obj.Prob.cas_opts = opts;
else
    error('Options passed to solver are of an unrecognized type.')
end
% selected solver-specific options
if ~isfield(obj.Prob.cas_opts,'options')
    obj.Prob.cas_opts.options = struct();
end

% available solvers -> qpsol or nlpsol
nlp_solvers = {'AmplInterface','blocksqp','bonmin','ipopt','knitro','snopt','worhp','scpgen','sqpmethod'};
qp_solvers  = {'cplex','gurobi','ooqp','qpoases','sqic','nlpsol'};

% regular problem
prob    = struct('f', obj.Prob.cost, 'x', obj.Prob.x_, 'g', [A;Adyn]*obj.Prob.x_,'p',obj.Prob.p_);
if contains(obj.Prob.cas_opts.solver,qp_solvers)
    obj.Prob.QP1 = casadi.qpsol( 'solver', obj.Prob.cas_opts.solver,prob,obj.Prob.cas_opts.options);
elseif contains(obj.Prob.cas_opts.solver,nlp_solvers)
    obj.Prob.QP1 = casadi.nlpsol('solver', obj.Prob.cas_opts.solver,prob,obj.Prob.cas_opts.options);
else
    error('Unrecognized solver specified');
end

% make get functions
obj.Prob.get_lba = casadi.Function('get_lba',{obj.Prob.p_},{[lba;ulba_dyn]});
obj.Prob.get_uba = casadi.Function('get_uba',{obj.Prob.p_},{[uba;ulba_dyn]});
get_lbx = casadi.Function('get_lbx',{obj.Prob.p_},{lbx});
obj.Prob.get_lbx = @(p) max(full(get_lbx(p)),[],2);
get_ubx = casadi.Function('get_ubx',{obj.Prob.p_},{ubx});
obj.Prob.get_ubx = @(p) min(full(get_ubx(p)),[],2);

% initializing result structure
zero_g = zeros(size(A,1)+size(Adyn,1),1);
obj.Prob.res = struct;
obj.Prob.res.x     = zero_x;
obj.Prob.res.lam_x = zero_x;
obj.Prob.res.g     = zero_g;
obj.Prob.res.lam_g = zero_g;
obj.Prob.get_x0     = @() obj.Prob.res.x(1:size(zero_x,1),1);    % indices specified to accomodate use after backup solver
obj.Prob.get_lam_x0 = @() obj.Prob.res.lam_x(1:size(zero_x,1),1);
obj.Prob.get_lam_a0 = @() obj.Prob.res.lam_g(1:size(zero_g,1),1);

obj.Prob.p2res = @(p) obj.Prob.QP1('p',p,...
    'x0',obj.Prob.get_x0(),...
    'lam_x0',obj.Prob.get_lam_x0(),'lam_g0',obj.Prob.get_lam_a0(),...
    'lbg',obj.Prob.get_lba(p),'ubg',obj.Prob.get_uba(p),...
    'lbx',obj.Prob.get_lbx(p),'ubx',obj.Prob.get_ubx(p));

obj.Prob.res2ufyf = @(res) deal(full(obj.Prob.x2uf(res.x)),full(obj.Prob.x2yf(res.x)));

%% create backup solver for when QP is infeasible
obj.Prob.backup = struct;
nxLarge = size(obj.Prob.x_,1); % number of optimization variables
na = size(A,1); % current number of constraints of form: lba <= Ax <= uba

% ---------- turning constraints soft (except for dynamics) ---------------
% relaxing constraints on x of form: lbx <= x <= ubx
ubx_zero = full(obj.Prob.get_ubx(zero_p)); % to see where inf
lbx_zero = full(obj.Prob.get_lbx(zero_p)); % to see where -inf 
mask_ulbx = ~isinf(lbx_zero) | ~isinf(ubx_zero);
np = sum(mask_ulbx);             % # of to be relaxed constraints
idx_relaxed = find(mask_ulbx); % row # of to be relaxed constraints

% combining all relaxed constraints of form: lba <= Ax <= uba
% -> new A matrix
P = speye(nxLarge,nxLarge); P = P(mask_ulbx,:);
npa= na+np;
A = [[P;A] speye(npa,npa)];
% -> new lba vector: [max(lbx(idx_relaxed,:),[],2),lba,ulba_dyn]
lbx_relaxed = lbx(idx_relaxed,:); % lbx has multiple columns!
get_lbx_relaxed = casadi.Function('get_lbx_relax',{obj.Prob.p_},{lbx_relaxed});
get_lbx_relaxed2= @(p) max(full(get_lbx_relaxed(p)),[],2);
get_lba_relaxed_part1 = casadi.Function('get_lba_relax_p1',{obj.Prob.p_},{[lba;ulba_dyn]});
obj.Prob.backup.get_lba = @(p) [get_lbx_relaxed2(p); full(get_lba_relaxed_part1(p))];
% -> new uba vector: [min(ubx(idx_relaxed,:),[],2),uba,ulba_dyn]
ubx_relaxed = ubx(idx_relaxed,:); % ubx has multiple columns!
get_ubx_relaxed = casadi.Function('get_ubx_relax',{obj.Prob.p_},{ubx_relaxed});
get_ubx_relaxed2= @(p) min(full(get_ubx_relaxed(p)),[],2);
get_uba_relaxed_part1 = casadi.Function('get_uba_relax_p1',{obj.Prob.p_},{[uba;ulba_dyn]});
obj.Prob.backup.get_uba = @(p) [get_ubx_relaxed2(p); full(get_uba_relaxed_part1(p))];

% enlarge optimization variable with slack variables to relax constraints
obj.Prob.backup.sigma_ = casadi.SX.sym('sigma', npa, 1);
obj.Prob.backup.x_     = vertcat(obj.Prob.x_,obj.Prob.backup.sigma_);

% enlarge matrix for dynamics constraints
if obj.options.ExplicitPredictor
    Adyn = [Adyn,sparse(obj.ny*obj.f,npa)];
else
    Adyn = [Adyn,sparse(obj.nGcols*(obj.nu+obj.ny)*obj.pfid,npa)];
end

% alter cost: heavily penalize constraint violations
obj.Prob.backup.cost = obj.Prob.cost+1e15*obj.Prob.backup.sigma_.'*obj.Prob.backup.sigma_;

% set up relaxed optimization problem
prob2    = struct('f', obj.Prob.backup.cost, 'x', obj.Prob.backup.x_, 'g', [A;Adyn]*obj.Prob.backup.x_,'p',obj.Prob.p_);
if contains(obj.Prob.cas_opts.solver,qp_solvers)
    obj.Prob.backup.QP2 = casadi.qpsol( 'solver',obj.Prob.cas_opts.solver,prob2,obj.Prob.cas_opts.options);
elseif contains(obj.Prob.cas_opts.solver,nlp_solvers)
    obj.Prob.backup.QP2 = casadi.nlpsol('solver',obj.Prob.cas_opts.solver,prob2,obj.Prob.cas_opts.options);
else
    error('Unrecognized solver specified');
end

% create simple function to get result from parameters
obj.Prob.backup.p2res = @(p) obj.Prob.backup.QP2('p',p,...
    'lbg',obj.Prob.backup.get_lba(p),'ubg',obj.Prob.backup.get_uba(p));
end

%% Helper functions
% =================== parsing user-defined constraints ====================
function usr_con = parse_usr_con(obj, usr_con,str1,str2)
% parses user-defined constraints:
% -> removes empty fields
% -> checks dimensions & possible classes

% possible fields for user-defined constraints:
all_str = [str1(:)',str2(:)',{'expr'}];

% get field names
fns = fields(usr_con);

% throw error if there are unrecognized constraint fields
fns_unrec = setdiff(fns,all_str);
if ~isempty(fns_unrec)
    % create string with unrecognized fields
    str_unrec = cell2mat( cellfun(@(x) [x,', '],fns_unrec,UniformOutput=false) );
    str_unrec = str_unrec(1:end-2);
    % create string with allowable fields
    str_rec = cell2mat( cellfun(@(x) [x,', '],all_str,UniformOutput=false) );
    str_rec = str_rec(1:end-2);
    error('Unrecognized constraint field(s): %s\nAllowed constraint fields are: %s',str_unrec ,str_rec);
end

% ---------------- for each of the entries in str1, st2: ------------------
% -> delete empty fields
% -> check class/dimensions

for k = 1:length(fns)
    fn = fns{k};          % field name
    value = usr_con.(fn); % field value 
    % -------------------- delete empty fields ------------------------
    if isempty(value)
        usr_con = rmfield(usr_con,fn);  % remove empty field
        % update field names at the end to prevent skipping others in loop
        continue;
    elseif strcmp(fn,'expr')
        continue;
    end
    
    % ------------------------ str1 & str2 ------------------------
    % => fields in str1 & str2 are of an allowed class
    % => fields in str1 & str2 are of the correct dimensions
    switch fn
        case {'uf','u0','u_min','u_max','du_max'}
            uy_dim = obj.nu;
        case {'yf','rf','y0','y_min','y_max','dy_max'}
            uy_dim = obj.ny;
    end
    switch fn
        case {'uf','yf','rf'}
            pos_classes = 'casadi.SX';
            pos_size    = [uy_dim, obj.f];
        case {'u0','y0'}
            pos_classes = 'casadi.SX';
            pos_size    = [uy_dim, 1];
        case {'u_min','u_max','y_min','y_max','du_max','dy_max'}
            pos_classes = 'double';
            pos_size    = [uy_dim, 1];
    end

    try
        validateattributes(value,{pos_classes},{'ndims',2,'size',pos_size})
    catch Error
        disp(append('For the field recognized as ',fn,' the following eror was encountered:'))
        error(Error.message)
    end

end
fns = fields(usr_con); % update field names (some may be removed)

% --------------------------- flag for str1 -------------------------------
if isempty(intersect(fns,str1))
   str1_flag = 'none';
else
   str1_flag = 'SX'; % casadi.SX is only other option
end

% ------------------------ expr handling and flag -------------------------
% possible flags: none, SX
if isfield(usr_con,'expr')
    if isa(usr_con.expr,'cell')
        if isvector(usr_con.expr)
            con_LHS  = usr_con.expr(1:2:end); con_LHS = con_LHS(:);
            con_gleq = usr_con.expr(2:2:end); con_gleq= con_gleq(:);
        elseif ismatrix(usr_con.expr)
            con_LHS  = usr_con.expr(:,1);
            con_gleq = usr_con.expr(:,2);
        else
            error('Incorrect cell array structure specified for the expr field.')
        end
        if length(con_LHS) ~= length(con_gleq)
            error('Incorrect cell array structure specified for the expr field.')
        elseif ~all(cellfun(@(x) isa(x,'char'),con_gleq))
            error('Incorrect cell array structure specified for the expr field.')
        end
        if all(cellfun(@(x) isa(x,'casadi.SX'),con_LHS))
            expr_flag = 'SX';
        else
            error('Left hand side of constraints may only contain expressions of casadi.SX class.')
        end
    else
        error('Constraint expression(s) must be supplied in a cell array.')
    end
else
    expr_flag = 'none'; % -> no 'expr' field
end

% -------------------------------- consistency between str1 & expr ------------------------
if strcmp(expr_flag,'SX') && strcmp(str1_flag,'none')
    % create string of all variable names in str1
    str1_str = cell2mat( cellfun(@(x) [x,', '],str1,UniformOutput=false) );
    str1_str = str1_str(1:end-2);
    error('Provide defined variables in expressions of expr field of constraints. Choose from: %s',str1_str);
elseif ~strcmp(expr_flag,str1_flag)
    error('Inconsistent use of specified parameter/variable type and expression.');
end
end