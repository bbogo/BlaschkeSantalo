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
