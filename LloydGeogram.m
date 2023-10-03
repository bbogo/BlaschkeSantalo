function str = LloydGeogram(samples,bbox,lb,ub,niter,pl)

if nargin<5
niter = 100;
end

ncells = length(samples);
fprintf("Number of cells: %d\n",ncells);

opt.ncells = ncells;
opt.bbox=bbox;
opt.dim = size(samples,1);
opt.plot=0;
	fun = @(x)BSMatFun(x,opt);
tol = 1e-4;

frac = 1;

Nfeval = 0;

for i=1:niter
	% compute Voronoi diagram and Centroids
	coords = zeros(ncells,2);
    for j=1:ncells
        coords(j,:) = BSMatFun(samples(:,j),opt);
    end    
    %if mod(i,50)==0
   % 	samples = BSMatReCenterAll(samples,-1,1);
   % end
    
    if nargin>=6
		[Areas, Centroids, Moments] = VoronoiGeogram(coords,bbox,1);
		%drawnow
		%pause
	else
		[Areas, Centroids, Moments] = VoronoiGeogram(coords,bbox);
	end
	%pause
	% move samples towards centroids
	mat = coords-Centroids;
	dists = sqrt(sum(mat.^2,2));
	[sd,ind] = sort(dists,'descend');
	val = sum(dists);

	newsam = samples;
	ncomp = 0;
	
	errs = zeros(ncells,1);
	
	omega = 1.5;

	for j=1:(frac*size(Centroids,1))
		%if norm(coords(j,:)-Centroids(j,:))>tol
		    ncomp = ncomp+1;
		    ap = BSMatFun(samples(:,ind(j)),opt);
			
			
			[newsam(:,ind(j)),out] = MoveCenter4(samples(:,ind(j)),Centroids(ind(j),:),@(x)BSMatFun(x,opt),@(x) BSMatJac(x,opt),lb,ub);
			
			Nfeval = Nfeval+out.funcCount;
			
			%newsam(:,ind(j)) = MoveCenter2(samples(:,ind(j)),ap+omega*(Centroids(ind(j),:)-ap),@(x)BSMatFun(x,opt),@(x) BSMatJac(x,opt),lb,ub);
		
		%else
		    %pause
	%		newsam(:,j) = samples(:,j);
	%	end
		errs(j) = norm(fun(newsam(:,ind(j)))-fun(samples(:,ind(j))));
	end
	dif = max(errs);
	fprintf('Iter %4d | Obj %.6f | Nopt %4d | VarObj %.6f | Diff %.6f | Nfeval %5d \n',i,val,ncomp,sum(Moments),dif,Nfeval);
	

	if dif<tol
		break
	end
	samples = newsam;
end


str.samples = samples;
str.bbox  = bbox;



function res = MoveCenter(x,c,fun,jac,lb,ub)

	%options.GradConstr = 'on';
	options.GradObj = 'on';
	options.Display = 'iter';
	options.MaxIterations = 200;
	options.algorithm = 'active-set';
	%options.Hessian = {'lbfgs',10};

	options.Display = 'none';
	options.ConstraintTolerance = 1e-6;
	options.StepTolerance = 1e-6;
	options.OptimalityTolerance = 1e-7;
	%options.CheckGradients='on';
	%options.UseParallel = true;
	%options.SubproblemAlgorithm = 'cg';


	vv  = fun(x);
	opt.dim = length(x);
	opt.vv  = vv;
	opt.fun = fun;
	opt.jac = jac;
	init = x;
	Lb = lb*ones(size(init));
	Ub = ub*ones(size(init));
	
	[res,fval,exflag,lambda,grad,hess] = fmincon(@(x) Obj(x,c,fun,jac),init,[],[],[],[],Lb,Ub,[],options);
    if exflag<1
    	fprintf("Exit flag %d\n",exflag);
    	pt = fun(res);
    	fprintf("Point %.6f %.6f\n",pt(1),pt(2));
    end
    
    %pause
    
function res = MoveCenter2(x,c,fun,jac,lb,ub)

	% knitro
	options = knitro_options('algorithm', 3, 'outlev', 0 , 'gradopt', 1, ...
                         'hessopt',6,'maxit', 1000, 'strat_warm_start',1,'honorbnds',1, ...
                         'opttol', 1e-8, ...
                         'lmsize' , 10,...
                         'presolve_initpt',0);
	% Specify an initial point.
	x0 = [-2; 1];

	% Specify some model characteristics.
	A = []; b = [];
	Aeq = []; beq = [];
	
	% Call Knitro to solve the optimization model.
	

	vv  = fun(x);
	opt.dim = length(x);
	opt.vv  = vv;
	opt.fun = fun;
	opt.jac = jac;
	init = x;
	Lb = lb*ones(size(init));
	Ub = ub*ones(size(init));
	
	[res,fval,exitflag,output,lambda,grad] = ...
    knitro_nlp(@(x) Obj(x,c,fun,jac),init,A,b,Aeq,beq,Lb,Ub,[],[],options);

    if exitflag<0
    	fprintf("Exit flag %d\n",exitflag);
    	pt = fun(res);
    	fprintf("Point %.6f %.6f\n",pt(1),pt(2));
    end
    
    %pause
function [res,output] = MoveCenter4(x,c,fun,jac,lb,ub)

	% knitro
	c = c(:);
	
	% Specify an initial point.
	x0 = [-2; 1];

	% Specify some model characteristics.
	A = []; b = [];
	Aeq = []; beq = [];
	
	% Call Knitro to solve the optimization model.
	

	vv  = fun(x);
	opt.dim = length(x);
	opt.vv  = vv;
	opt.fun = fun;
	opt.jac = jac;
	init = x;
	Lb = lb*ones(size(init));
	Ub = ub*ones(size(init));
	

	%options.Algorithm = 'levenberg-marquardt';
	options = optimoptions('lsqnonlin');
    %options.Algorithm = 'levenberg-marquardt';
	options.SpecifyObjectiveGradient = true;
	%[res,fval,exitflag,output,lambda,grad] = ...
    %knitro_nlp(@(x) Obj(x,c,fun,jac),init,A,b,Aeq,beq,Lb,Ub,[],[],options);
   	options.Display='none';

    fff = @(x) ObjVec(x,c,fun,jac);
    [res,resnorm,resid,exitfl,output] = lsqnonlin(fff,init,Lb,Ub,options);
	
	
	
	resnorm = 1;
	if and(1==0,resnorm>1e-3)
		fprintf("Value = %.6f\n",resnorm);
		hold on
		plot(c(1),c(2),'.g');
		res
		drawnow
		hold off
	end
    %if exitflag<0
   % 	fprintf("Exit flag %d\n",exitflag);
   % 	pt = fun(res);
   % 	fprintf("Point %.6f %.6f\n",pt(1),pt(2));
   % end
    
    %pause
	
function [val,grad] = Obj(x,c,fun,jac)

pt = fun(x);
jj = jac(x);


val = sum((pt-c).^2);
grad = 2*jj.'*(pt(:)-c(:));


function [val,grad] = ObjVec(x,c,fun,jac)
epst = 1e-3;

pt = fun(x);
jj = jac(x);

%pt
%c


val = pt(:)-c(:);
grad = jj;


val = [val; epst*x(:).^2];
grad = [grad;
        epst*spdiags(x(:),0,length(x(:)),length(x(:)))
         ];


