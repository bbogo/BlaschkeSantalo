function res = BSrefine(str,nadd,rad)
str
try
samples = str.samples;
dim     = str.dim;
catch
samples = str.cells;
dim     = size(samples,1);
end

opt.dim = dim;
fun     = @(x)APWfunc(x);
jac     = @(x)APWjac(x);

%res = AddPoints(samples,nadd,fun,jac,1e-4,20,rad);
%res = AddPoints2(samples,nadd,fun,jac,1e-4,20,rad);
res = AddPoints3(samples,nadd,fun,jac,1e-4,20);
%max(res(:))
%min(res(:))
