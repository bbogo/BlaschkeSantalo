function newpts = AddPoints3(samples,nadd,fun,jac,lb,ub)

pl = 1;

ncells = size(samples,2);
dim    = size(samples,1);


newpts = samples;

thetas = linspace(0,2*pi,nadd+1);
thetas = thetas(1:end-1); 

coords = zeros(ncells,2);
if pl==1
figure(1)
clf
hold on
end
% Compute coordinates
for i=1:ncells
        coords(i,:) = fun(samples(:,i));
end
if(pl==1)
plot(coords(:,1),coords(:,2),'.b');
hold off
end

% search nearest neighbors
Idx = knnsearch(coords,coords,'K',nadd,'IncludeTies',false);
dt = delaunayTriangulation(coords(:,1),coords(:,2));
tri = dt.ConnectivityList;
edges = [tri(:,1) tri(:,2); 
         tri(:,2) tri(:,3); 
         tri(:,3) tri(:,1)];
edges    = sort(edges,2);
edges    = unique(edges,'rows');


ledges = sqrt(sum((coords(edges(:,1),:)-coords(edges(:,2),:)).^2,2));
meanl   = mean(ledges);


% add midpoints
for i=1:size(edges,1)
	
	    if and(ledges(i)<1.75*meanl,ledges(i)>0.25*meanl)
			cmid = 0.5*(coords(edges(i,1),:)+coords(edges(i,2),:));
			res  = MoveCenter4(samples(:,edges(i,1)),cmid,fun,jac,lb,ub);
			newpts = [newpts res(:)];
			if(pl==1)
				hold on
				pt = fun(res);
				plot(coords(edges(i,:),1),coords(edges(i,:),2),'b');
				plot(pt(1),pt(2),'.r');
				
				hold off
			end
		end
	
end

%axis equal

%size(newpts)

function res = MoveCenter(x,c,fun,jac,lb,ub)

	%options.GradConstr = 'on';
	options.GradObj = 'on';
	options.Display = 'iter';
	options.MaxIterations = 50;
	options.algorithm = 'interior-point';
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
	c = c(:);
	options = knitro_options('algorithm', 3, 'outlev', 0 , 'gradopt', 1, ...
                         'hessopt',0,'maxit', 100, 'strat_warm_start',1,'honorbnds',1, ...
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

    %if exitflag<0
   % 	fprintf("Exit flag %d\n",exitflag);
   % 	pt = fun(res);
   % 	fprintf("Point %.6f %.6f\n",pt(1),pt(2));
   % end
    
    %pause
    
function res = MoveCenter3(x,c,fun,jac,lb,ub)

	% knitro
	c = c(:);
	options = knitro_options('outlev', 0 , 'gradopt', 1, ...
                         'maxit', 50, 'strat_warm_start',1,'honorbnds',1, ...
                         'opttol', 1e-8, ...
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
	
	%[res,fval,exitflag,output,lambda,grad] = ...
    %knitro_nlp(@(x) Obj(x,c,fun,jac),init,A,b,Aeq,beq,Lb,Ub,[],[],options);
    [res,resnorm,residual,exitflag,output,lambda,jacobian] = ...
    knitro_nlnlsq(@(x) ObjVec(x,c,fun,jac),init,Lb,Ub,[],options);
	
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
function res = MoveCenter4(x,c,fun,jac,lb,ub)

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
    res = lsqnonlin(fff,init,Lb,Ub,options);
	
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

try
  [pt,jj] = fun(x);
catch
  pt = fun(x);
  jj = jac(x);
end

val = sum((pt(:)-c(:)).^2);
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

