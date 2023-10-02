function res = BSReCenter(sample,fun,jac,lb,ub)

% minimize |x|^2
% such that f(x) = f(sample)
% i.e. go on the level set, preserving the constraint

meth=1;

switch meth
case 1

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
case 2
	options = knitro_options('algorithm', 3, 'outlev', 0 , 'gradopt', 1,'gradconstr',1, ...
                         'hessopt',0,'maxit', 100, 'strat_warm_start',1,'honorbnds',1, ...
                         'opttol', 1e-8, ...
                         'lmsize' , 10,...
                         'presolve_initpt',0);
                         
    vv  = fun(sample);
	opt.dim = length(sample);
	opt.vv  = vv;
	opt.fun = fun;
	opt.jac = jac;
	% Specify some model characteristics.
	A = []; b = [];
	Aeq = []; beq = [];
	
	init = sample(:);
	Lb = lb*ones(size(init));
	Ub = ub*ones(size(init));
	

	
	[res,fval,exitflag,output,lambda,grad] = ...
    knitro_nlp(@(x) Obj(x),init,A,b,Aeq,beq,Lb,Ub,[],@(x) nonlcon(x,opt),options);
    res
end	


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


