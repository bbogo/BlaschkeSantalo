function [Areas,Centroids,Moments,mind] =  VoronoiGeogram(pts,bbox,plotting)
% input Voronoi points
% points of bounding box
% triangulation of bounding box



N = size(pts,1);

if nargout>3
	% compute minimal distances, nearest neighbors
	Idx = knnsearch(pts,pts,'K',2);
	mind = zeros(size(pts,1));
	nn   = pts(Idx(:,2),:);
	%size(
	mind = sqrt(sum((pts-nn).^2,2));
end

% save bounding box and points for Geogram
%meshio.write_points_cells("box.mesh",points,tri)

writeMesh(bbox);


dlmwrite('points.xyz',pts,'delimiter',' ','precision','%.16e');

%tic
% run Geogram
geos = '~/geogram/build/Linux64-gcc-dynamic-Release/bin/compute_RVD ';
str = [geos,' ./box.mesh ./points.xyz ./voronoioutput.obj epsilon=0 compression_level=0 use_doubles=true algo:predicates=exact integration_smplx=false  cell_borders=true volumetric=false '];
[status,res]  = unix(str);

%toc

[s,r]=system("sed -n -e 's/^f//p' voronoioutput.obj > faces.txt");
[s,r]=system("sed -n -e 's/^v//p' voronoioutput.obj > vertices.txt");
[s,r]=system("sed -n -e 's/^# attrs f //p' voronoioutput.obj > attributes.txt");  

% load vertices
vert = dlmread("vertices.txt");
vert = vert(:,1:2);

% load attributes
attrs = dlmread("attributes.txt");
attrs(:,2) = attrs(:,2)+1;


% load faces
% result is array padded with zeros... good!
faces = dlmread("faces.txt");

if nargin>2
	try 
    	col = plotting.col;
    catch
    	col = rand(N,3);
    end
    figure(1)
    set(gcf,'DefaultAxesFontSize',18,'DefaultAxesTickLabelInterpreter','latex')
    clf
    hold on
end

Areas     = zeros(N,1);
Meanxy    = zeros(N,2);
Centroids = zeros(N,2);
Jx        = zeros(N,1);
Jy        = zeros(N,1);
Moments   = zeros(N,1);


for i=1:length(faces)
	face = flip(faces(i,find(faces(i,:)>0)));
	
	Ip = attrs(i,2);

	xx = vert(face,1);
	yy = vert(face,2);
	
	xp = circshift(xx,-1);
	yp = circshift(yy,-1);
	
	if nargin>2
	%0.6+0.4*col(Ip,:)
	fill([xx;xx(1)],[yy;yy(1)],0.6+0.4*col(Ip,:),'EdgeColor',[0.5,0.5,0.5],'LineWidth',1.5)
	%plot([xx;xx(1)],[yy;yy(1)],'b')
	end
	Areas(Ip) = Areas(Ip)+0.5*(dot(xx,yp)-dot(xp,yy));
	% Centroids (do not divide by the area now)
    Meanxy(Ip,1) = Meanxy(Ip,1)+1/6*dot((xx+xp),(xx.*yp-xp.*yy));
    Meanxy(Ip,2) = Meanxy(Ip,2)+1/6*dot((yy+yp),(xx.*yp-xp.*yy));


    % Moments
    Jx(Ip) = Jx(Ip)+1/12*dot((yy.^2+yy.*yp+yp.^2),(xx.*yp-xp.*yy));
    Jy(Ip) = Jy(Ip)+1/12*dot((xx.^2+xx.*xp+xp.^2),(xx.*yp-xp.*yy));
	
	
end


%Areas
Centroids = Meanxy;
Centroids = Meanxy./repmat(Areas(:),[1,2]);
Cplot     = Centroids;
Moments   = Jx+Jy-2*pts(:,1).*Meanxy(:,1)-2*pts(:,2).*Meanxy(:,2)+Areas.*(pts(:,1).^2+pts(:,2).^2);

if nargin>2

	msize = 10;
	if N>200
		msize = 8;
	end
	if N>500
		msize = 5;
	end
	if N>2000
		msize = 5;
	end



plot(pts(:,1),pts(:,2),'.','color',[1,0,0],'MarkerSize',msize);%,'MarkerEdgeColor','k',...
    %'MarkerFaceColor','g','MarkerSize',2);
    
plot(Cplot(:,1),Cplot(:,2),'.b','MarkerSize',msize);

Ax = gca;

switch bbox(1,2)
	case -2.5
	Ax.YTick = [-2,0,2];
	Ax.XTick = [-2,0,2];
	case -5
	Ax.XTick = [-3,0,3];
	Ax.YTick = [-4,0,4];
	case -20
	Ax.XTick = [-4,0,4];
	Ax.YTick = [-20,-10,0,10,20];
end

Ax.XGrid = 'on';
Ax.YGrid = 'on';
Ax.Layer = 'top';
Ax.GridLineStyle = '-';
Ax.GridAlpha = 0.5;

else
	plot(pts(:,1),pts(:,2),'sr','MarkerSize',15,'MarkerEdgeColor','red',...
    'MarkerFaceColor',[1 .6 .6]);
	plot(Cplot(:,1),Cplot(:,2),'.b','MarkerSize',15);

end

%dt = delaunayTriangulation(pts(:,1),pts(:,2));
%triplot(dt);


hold off
axis equal 
axis tight
%axis off









