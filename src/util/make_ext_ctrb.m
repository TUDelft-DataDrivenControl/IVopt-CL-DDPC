function Ks = make_ext_ctrb(A,B,s,options)
% make extended (reversed) controllability matrix
arguments
    A double
    B double
    s double {isscalar}
    options.rev logical = false; % -> set to true if reversed
end
Ks = cell(1,s);

if ~options.rev
    Ks{1} = B;
    for k = 2:s
        Ks{k} = A * Ks{k-1};
    end
else
    Ks{s} = B;
    for k = s:-1:2
         Ks{k-1} = A  *  Ks{k};
    end
end
Ks = cell2mat(Ks);
end