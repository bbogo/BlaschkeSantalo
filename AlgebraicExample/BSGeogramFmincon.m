function str = BSGeogram(samples,bbox)

	% compute restricted Centroidal Voronoi Diagram
	% for the image 
	
	% Plotting (or not) during optimization
	% Choosing 1 slows down the code
	opt.plot=0;
	


    rad = 0.5;
    cen = [2,1];

    % turn samples into points in R^d

    ncells = size(samples,2);
    fprintf("Number of cells: %d\n",ncells);
    
    %for it in range(0,nitLloyd):
    %    ps, bnd, areas, vor, ridge_points = PolyVoronoi(coords,poly,debug=debug,fact=fact)
    %    for i in range(0,len(coords)):
    %        coords[i] = ps[i].centroid.coords
    % use scipy minimize
    
    
    
    
    init = samples(:);
    d    = rand(size(init));
    %d(:) = 0;
    %d(5) = 1;
    epst = 1e-4;
    
    opt.ncells = ncells;
    
    opt.col = rand(ncells,3); % random colors for plot
    
    opt.bbox=bbox;
    opt.dim = size(samples,1);
    [val,grad] = ObjFunc(init,opt);

    v1 = ObjFunc(init+epst*d,opt);
    v2 = ObjFunc(init-epst*d,opt);
    
    fprintf("Analytic gradient: %e\n",dot(grad,d))
    fprintf("Finite diff:       %e\n",(v1-v2)/2/epst)
    
	lb = -1*ones(size(init));
	ub =    ones(size(init));
	
	options = optimoptions('fmincon','Algorithm','active-set');
	options.SpecifyObjectiveGradient = true;
	options.Display = 'iter';
	
	% Modify this for more precise results
	options.MaxFunctionEvaluations = 200;
	options.MaxIterations = 100;
	options.HessianApproximation = {'lbfgs',15};
	options.FunctionTolerance = 1e-8;
	options.OptimalityTolerance = 1e-6;


	res = fmincon(@(x) ObjFunc(x,opt),init,[],[],[],[],lb,ub,[],options);

	opt.plot=1;
	[val, gradt,mind]=ObjFunc(res,opt);

	str.samples = reshape(res,[],ncells);
	str.bbox = bbox;
	str.ncells = ncells;
	str.dim  = opt.dim;
	str.mind = mind;

    
    
function [val, gradt,mind] = ObjFunc(x,opt)
	
	

    ncells = opt.ncells;    
    bbox   = opt.bbox;
    mat  = reshape(x,opt.dim,[]);
        
    coords = zeros(ncells,2);
    for i=1:ncells
        coords(i,:) = BSMatFun(mat(:,i),opt);
    end    
        
        if opt.plot==0
        	[Areas, Centroids, Moments] = VoronoiGeogram(coords,bbox);
		else
			[Areas, Centroids, Moments,mind] = VoronoiGeogram(coords,bbox,opt);
			drawnow
		end
        
        val = sum(Moments);
        
        gradt = zeros(size(mat,1),ncells);

        for i= 1:ncells
            actpt = mat(:,i);
            actJac = BSMatJac(actpt,opt);
            gradt(:,i) = 2*Areas(i)*actJac.'*(coords(i,:)-Centroids(i,:)).';
        end
     
        gradt = gradt(:);
        
        
    
    


    
