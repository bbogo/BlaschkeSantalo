function res = BSMatReCenter(samples,lb,ub)

% minimize |x|^2
% such that f(x) = f(sample)
% i.e. go on the level set, preserving the constraint



%samples(:,1)

opt.dim = size(samples,1);

res = zeros(size(samples));

for i=1:size(samples,2)

	res(:,i) = APWReCenter(samples(:,i),@(x)APWfunc(x),@(x) APWjac(x),lb,ub);
end	


