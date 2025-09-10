classdef IV_4_DDPC < dynamicprops
    %IV_4_DDPC Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        p (1,1) double {mustBeInteger}
        f (1,1) double {mustBeInteger}
        nu (1,1) double {mustBeInteger}
        ny (1,1) double {mustBeInteger}
        N (1,1) double {mustBeInteger}
        Nbar (1,1) double {mustBeInteger}
        IV_names cell
        IV_descr cell
    end
    properties (Hidden)
        u double {mustBeNumeric,mustBeMatrix}
        y double {mustBeNumeric,mustBeMatrix}
    end
    properties (Dependent,Hidden) % noisy data
        Up
        Uf
        Yp
        Yf
        Wp % [Up; Yp]
    end
    properties (Dependent) % IVs
        iv1
    end
    
    methods
        function obj = IV_4_DDPC(u,y,p,f) % class constructor

            % set properties
            [obj.u,obj.y,obj.p,obj.f] = deal(u,y,p,f); % set properties
            [obj.nu,obj.Nbar] = size(u);
            obj.ny = size(y,1);
            obj.N = obj.Nbar - (p+f) + 1;

            % validate size of data
            validateattributes(y,{'double'},{'size',[NaN,obj.Nbar]});

            % create cell array with IV names
            obj.IV_names{1} = 'iv1';
            
            % add description for open-loop IV
            obj.IV_descr{1} = 'open-loop IV';
        end
        
        %% Get and Set methods (excl. IV matrices)
        function value = get.Up(obj)
            value = obj.make_Hankel(obj.u,obj.p);
            value = value(:,1:obj.N);
        end
        function value = get.Yp(obj)
            value = obj.make_Hankel(obj.y,obj.p);
            value = value(:,1:obj.N);
        end
        function value = get.Uf(obj)
            value = obj.make_Hankel(obj.u(:,obj.p+1:end),obj.f);
        end
        function value = get.Yf(obj)
            value = obj.make_Hankel(obj.y(:,obj.p+1:end),obj.f);
        end
        function value = get.Wp(obj)
            value = [obj.Up;obj.Yp];
        end

        %% Get and Set methods (IV matrices)
        % iv1: open-loop IV
        function Z = get.iv1(obj)
            Z = [obj.Wp;obj.Uf];
        end

        % method to add other IVs
        function add_IV(obj,IVname,Z0,opt)
            arguments
                obj
                IVname char
                Z0 double {mustBeMatrix}
                opt.TSLS logical = false;
                opt.descr char = '';
            end
            
            % check dimensions of Z0
            validateattributes(Z0,{'double'},{'size',[NaN,obj.N]});

            IVname_ = [IVname,'_']; % hidden IV property name
            % check whether property exists
            if isprop(obj,IVname)
                error("Property %s already exists",IVname)
            elseif isprop(obj,IVname_)
                error("Hidden property %s already exists",IVname_)
            end

            % add properties
            P1=addprop(obj,IVname);  % visible property
            P2=addprop(obj,IVname_); % hidden property
            P2.Hidden = true;
            
            % perform 2SLS?
            if opt.TSLS
                Z0 = obj.TSLS(obj.Uf,Z0);
            end

            obj.(IVname_) = Z0; % set hidden property
            P1.GetMethod = @(obj) [obj.Wp; obj.(IVname_)];

            % add name and description
            obj.IV_names{1,length(obj.IV_names)+1} = IVname;
            obj.IV_descr{1,length(obj.IV_descr)+1} = opt.descr;
        end
        
        % method to get description from name
        function descr = get_descr(obj,name)
            idx = find(strcmp(obj.IV_names,name),1,'first');
            if isempty(idx)
                warning('The IV name %s does not exist.\nReturning an empty description.',name);
                descr = '';
            else
                descr = obj.IV_descr{idx};
            end
        end

    end
    %% Static methods
    methods (Static)
        function Znew = TSLS(W,Z0)
            [Q, ~] = qr(Z0', 'econ'); % economy QR of Z0'
            % Q*Q.' is orthogonal projector onto the rowspace of Z0
            Znew = W*(Q*Q');
        end

        function Hs = make_Hankel(data, s)
        % Efficiently constructs past-future Hankel matrices from input data
        % Inputs:
        %   data: (ndata x Nbar) time series data matrix
        %   s: number of block rows
        % Outputs:
        %   Hs: (ndata*s x (Nbar - s + 1)) Hankel matrix
        
        [ndata, Nbar] = size(data);
        num_cols = Nbar - s + 1;
        
        Hs = zeros(ndata * s, num_cols);
        
        for i = 1:s
            Hs((i-1)*ndata+1:i*ndata, :) = data(:, i:i+num_cols-1);
        end

        end
    end
end

