function Z = monomials_exponents(vartable,d)
% MONOMIALS --- Construct a matrix of exponents of monomials with 
%       prespecified degrees.


i = find(d<0);
if ~isempty(i)
    error('Negative degree.');
end;

Z = [];
for i = d
    Z = [Z; oldconstructZ_exponents(vartable,i)];
end;

end

% ==================================================================
function Z = oldconstructZ_exponents(vartable,d)
% Old constructZ

if isnumeric(vartable)
    n = vartable;
else
    n = length(vartable);
end;

ZZ = sparse(1,n);
for i = 1:n
    ss = size(ZZ,1);
    ZZ = sprepmat(ZZ,d+1,1);
    for j = 0:d
        ZZ(ss*j+1:ss*j+ss,i) = j;
    end;
    idx = find(sum(ZZ,2) <= d);   % Throw away invalid monomials
    ZZ = ZZ(idx,:);
end;
idx = find(sum(ZZ,2) == d);
Z = ZZ(idx,:);

end
