function res = BSMatReCenter(samples,lb,ub)

% minimize |x|^p
% such that f(x) = f(sample)
% i.e. go on the level set, preserving the constraint


%samples(:,1)

opt.dim = size(samples,1);
res = zeros(size(samples));

for i=1:size(samples,2)
	res(:,i) = BSReCenter(samples(:,i),@(x)BSMatFun(x,opt),@(x) BSMatJac(x,opt),lb,ub);
end	


