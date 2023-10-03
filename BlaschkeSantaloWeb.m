function final=  BlaschkeSantaloWeb(N,d,nref)

% Computation of Blaschke Santalo diagram for the algebraic example
% Prerequisite: https://github.com/BrunoLevy/geogram
% Voronoi Diagrams are generated in VoronoiGeogram.m
% Modify the path to the 'compute_RVD' routine, if needed. 
% 
% Inputs:
% N = number of initial samples
% d = dimension (2,3,4), higher costs more
% nref = number of refinement (0,1,2,3... Do not exaggerate)
%
% Examples of usage:
% BlaschkeSantaloWeb(100,2,1)
% BlaschkeSantaloWeb(20,3,3)
% BlaschkeSantaloWeb(50,4,2)



setenv('LD_LIBRARY_PATH');
if nargin<1
	N = 10;
end


% dimension for the matrix problem
dim = d*(d+1)/2;

% choose bounding box
switch d
case 2
	ax = 2.5;
	ay = 2.5;
	bbox = [-ax -ay; ax -ay; ax ay; -ax ay];
case 3
	ax = 5;
	ay = 5;
	bbox = [-ax -ay; ax -ay; ax ay; -ax ay];
case 4
    ax = 6;
    ay = 20;
    bbox = [-ax -ay; ax -ay; ax ay; -ax ay];
end

% Write Box in Mesh format for Geogram
writeMesh(bbox);

% Initialization (random)
samples = 2*rand(dim,N)-1;
final.init = samples;
final.bbox = bbox;
final.cells = {};
pos = 1;

% First run
str = BSGeogramFmincon(samples,bbox);

fprintf("First optimization result.\n");
pause

final.cells{pos} = str.samples;

% Number of points to add on refinement
nadd = 5;

for i=1:nref
	% Re-centering
	fprintf("Recentering number %d\n",i);
	tic
	cen = BSMatReCenterAll(str.samples,-1,1);
	toc
	str.samples = cen;
	% Refinement
	fprintf("Refinement number %d\n",i);
	tic
	res = BSrefine(str,nadd,0.25*str.mind);
	toc
	pause
	
	% Run Lloyd or not
	%str = LloydGeogram(res,bbox,-1,1,100,1);
	str = BSGeogramFmincon(res,bbox);
	fprintf("Optimization after refinement number %d\n",i);
	
	final.cells{i+pos} = str.samples;
	final.mind         = str.mind;
	final.dim          = str.dim;
end

save('currentBS.mat','final');


function writeMesh(bbox)

files = 'box.mesh';

tri = delaunay(bbox);

bbox = [bbox ones(size(bbox,1),1)];
tri  = [tri ones(size(tri,1),1)];
fid = fopen(files,'w');
fprintf(fid,'MeshVersionFormatted 2\n');
fprintf(fid,'Dimension 2\n');
fprintf(fid,'\n');
fprintf(fid,'Vertices\n');
fprintf(fid,'%d\n',size(bbox,1));
fclose(fid);
dlmwrite(files,bbox,'delimiter',' ','-append');


fid = fopen(files,'a');
fprintf(fid,'Triangles\n');
fprintf(fid,'%d\n',size(tri,1));
fclose(fid);
dlmwrite(files,tri,'delimiter',' ','-append');
fid = fopen(files,'a');
fprintf(fid,'End\n');
fclose(fid);
