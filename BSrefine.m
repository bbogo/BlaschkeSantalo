function res = BSrefine(str,nadd,rad)

try
samples = str.samples;
dim     = str.dim;
catch
samples = str.cells;
dim     = size(samples,1);
end

opt.dim = dim;
fun     = @(x)BSMatFun(x,opt);
jac     = @(x)BSMatJac(x,opt);

%res = AddPoints(samples,nadd,fun,jac,-1,1,rad);
%res = AddPoints2(samples,nadd,fun,jac,-1,1,rad);
res = AddPoints3(samples,nadd,fun,jac,-1,1);


%max(res(:))
%min(res(:))
