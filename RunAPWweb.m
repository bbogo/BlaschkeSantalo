function final=  RunAPWweb(N,dim,nref)

% Computation of Blaschke Santalo diagram for the geometric example
% Prerequisite: https://github.com/BrunoLevy/geogram
% Voronoi Diagrams are generated in VoronoiGeogram.m
% Modify the path to the 'compute_RVD' routine, if needed. 
% 
% Inputs:
% N = number of initial samples
% d = number of parameters (50 works fine)
% nref = number of refinements (up to 3 the costs are reasonable)
%
% Examples of usage:
% RunAPWweb(15,50,2)
%
% Don't start with many samples as the'll all concentrate around the origin.




setenv('LD_LIBRARY_PATH');
if nargin<1
	N = 10;
end

NiterLloyd = 100;

figure(1);

bbox = [-1 -1; 9 -1; 9 7; -1 7];

writeMesh(bbox);


	samples = rand(dim,N);
	final.init = samples;
	final.bbox = bbox;
	final.cells = {};
	pos = 1;

col = rand(size(samples,1),3);
opt.col = col;
size(col)

fprintf("Initial run for Lloyd Algorithm\n");
tic

% without plotting (faster)
%str = APWLloydGeogram(samples,bbox,0,1,NiterLloyd);
% with plotting (slower)
str = APWLloydGeogram(samples,bbox,0,1,NiterLloyd,opt);
toc
pause
samples = str.samples;
nadd = 5;
final.cells{1} = str.samples;
save('tempAPW.mat','final');
for n=1:nref
	fprintf("Recentering CVT\n");
	tic
	cen = APWReCenterAll(str.samples,0,1);str.samples = cen;
	toc
	fprintf("Refining CVT\n");
	tic
	res = APWrefine(str,nadd);
	toc
	fprintf("Running Lloyd: refinement %d \n",n);
	tic
	col = rand(size(res,1),3);
	opt.col = col;
	% without plotting (faster)
	%str = APWLloydGeogram(res,bbox,0,1,NiterLloyd);
	% with plotting (slower)
	str = APWLloydGeogram(res,bbox,0,1,NiterLloyd,opt);
	toc
	final.cells{n+1} = str.samples;
	save('tempAPW.mat','final');
	final.cells{n+1} = str.samples;
	save('tempAPW.mat','final');
end

save('currentAPW.mat','final');


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
