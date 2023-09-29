function res = BSReCenter(sample,fun,jac,lb,ub)

% minimize |x|^2
% such that f(x) = f(sample)
% i.e. go on the level set, preserving the constraint



	options.GradConstr = 'on';
	options.GradObj = 'on';
	options.Display = 'iter';
	options.MaxFunEvals = 100;
	options.algorithm = 'interior-point';
	options.Hessian = {'lbfgs',10};

	options.Display = 'none';
	options.ConstraintTolerance = 1e-6;
	options.StepTolerance = 1e-6;
	options.OptimalityTolerance = 1e-5;
	%options.UseParallel = true;
	%options.SubproblemAlgorithm = 'cg';


	vv  = fun(sample);
	opt.dim = length(sample);
	opt.vv  = vv;
	opt.fun = fun;
	opt.jac = jac;
	init = sample;
	Lb = lb*ones(size(init));
	Ub = ub*ones(size(init));
	
	res = fmincon(@ Obj,init,[],[],[],[],Lb,Ub,@(x) nonlcon(x,opt),options);
	
function [val,grad] = Obj(x)
p = 10;

valp = sum(x.^p);

val  = valp^(1/p);
grad = x.^(p-1)*valp^(1/p-1);	
	
function [c,ceq,gradc,gradceq]	= nonlcon(x,opt)

c = [];
gradc = [];
ceq     = opt.fun(x)-opt.vv;
gradceq = opt.jac(x).';


