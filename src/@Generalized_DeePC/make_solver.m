function make_solver(obj,usr_con)

%% parsing user-defined constraints
usr_con = parse_usr_con(obj, usr_con);

%% specify how to make parameters/variables
obj.make_var = @(dim1,dim2,name) casadi.SX.sym(name,dim1,dim2);
obj.make_par = @(dim1,dim2,name) casadi.SX.sym(name,dim1,dim2);

%% make variables/parameters up, yp, rf, uf, yf

if isfield(usr_con,'u0')
    obj.Prob.up_ = [obj.make_par(obj.nu,obj.p-1,'up_endmin1'),usr_con.u0];
else
    obj.Prob.up_ = obj.make_par(obj.nu,obj.p,'up');
end
if isfield(usr_con,'y0')
    obj.Prob.yp_ = [obj.make_par(obj.ny,obj.p-1,'yp_endmin1'),usr_con.y0];
else
    obj.Prob.yp_ = obj.make_par(obj.ny,obj.p,'yp');
end
if isfield(usr_con,'rf')
    obj.Prob.rf_ = usr_con.rf;
else
    obj.Prob.rf_ = obj.make_par(obj.ny,obj.f,'rf');
end
if isfield(usr_con,'uf')
    obj.Prob.uf_ = usr_con.uf;
else
    obj.Prob.uf_ = obj.make_var(obj.nu,obj.f,'uf');
end
if isfield(usr_con,'yf')
    obj.Prob.yf_ = usr_con.yf;
else
    obj.Prob.yf_ = obj.make_var(obj.ny,obj.f,'yf');
end

%% Construct solver
obj.make_CasADi_solver(usr_con)

end

%% Helper functions
% =================== parsing user-defined constraints ====================
function usr_con = parse_usr_con(obj, usr_con)
% parses user-defined constraints:
% -> removes empty fields
% -> checks dimensions & possible classes

% two types of constraint specifications possible using below fields
str1 = {'uf','yf','rf','u0','y0'}; %-> fns should contain expr field too if used
str2 = {'u_min','u_max','y_min','y_max','du_max','dy_max'};
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