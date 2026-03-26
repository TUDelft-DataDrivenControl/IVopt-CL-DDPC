function Gs = make_ext_obsv(A,C,s)
% make extended observability matrix
Gs = cell(s,1);
Gs{1} = C;
for k = 2:s
    Gs{k}  = Gs{k-1}*A;
end
Gs  = cell2mat(Gs);
end