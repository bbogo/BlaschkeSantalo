function [A,P,W,GradA,GradP,GradW] = BSAPWalt(x)

% Alternate discretization of convex sets
% using second order differences as parameters


% symmetrize, find values of support function

N = length(x)-1;
w = x(1);     % width
k = x(2:end); % second order differences

[ix,iy] = meshgrid(1:N,1:N);
AA = N+1-max(ix,iy);

% direct parametrization xx,yy...
xx = linspace(0,1,N+1);
yy = w+[AA*k; 0];

xx = xx(:);
yy = yy(:);

% add points [0,0] and [1,0]
% to create a quarter polygon
%xx = [0; xx(:); 1];
%yy = [0; yy(:); 0];


%clf
%plot(xx,yy,'.r');
%axis equal

%axis([0 1.1 0 max(yy)+0.1])

%pause
h = 1/N;
A = 0.5*h*sum(yy(1:end-1)+yy(2:end));
ll = sqrt((yy(1:end-1)-yy(2:end)).^2+h^2);
P = w+sum(ll);

xp = xx(2:end);
xm = xx(1:end-1);
yp = yy(2:end);
ym = yy(1:end-1);

W = 1/12*(3*w+w^3+dot(xm.^2+xm.*xp+xp.^2+ym.^2+ym.*yp+yp.^2,xp.*ym-xm.*yp));

% Gradient Area OK
GradA = zeros(size(x));
GradA(1) = 1;
GradA(2:end) = 0.5*h*(sum(AA(:,1:end-1),2)+sum(AA(:,2:end),2)+1);

% Gradient Perim OK
GradP = zeros(size(x));
GradP(1) = 1;
B = (ix>=iy);
GradP(2:end) = B*((yy(1:end-1)-yy(2:end))./ll);

% Gradient Moment!


GradW = zeros(size(x));
%dot(xm.^2+xm.*xp+xp.^2+ym.^2+ym.*yp+yp.^2,xp.*ym-xm.*yp)
GradW(1) = 1/12*(3+3*w^2+dot(3*ym+3*yp,xp.*ym-xm.*yp)+...
                        +dot(xm.^2+xm.*xp+xp.^2+ym.^2+ym.*yp+yp.^2,xp-xm));

partW1 = (2*ym+yp).*(xp.*ym-xm.*yp)+(xm.^2+xm.*xp+xp.^2+ym.^2+ym.*yp+yp.^2).*(xp);

partW2        = (2*yp+ym).*(xp.*ym-xm.*yp)+(xm.^2+xm.*xp+xp.^2+ym.^2+ym.*yp+yp.^2).*(-xm);
partW1(2:end) = partW1(2:end)+partW2(1:end-1);

GradW(2:end) = 1/12*AA.'*partW1;

A = 4*A;
P = 4*P;
W = 4*W;

GradA = 4*GradA;
GradP = 4*GradP;
GradW = 4*GradW;
