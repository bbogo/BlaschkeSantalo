function Jac = APWjac(x)
scalex = 100;

[A,P,W,GradA,GradP,GradW] = BSAPWalt(x);

im = zeros(2,1);
im(1) = scalex*A/P^2;
im(2) = A^2/W;

Jac = zeros(2,length(x));
Jac(1,:) = scalex/P^2*GradA+scalex*A*(-2)/P^3*GradP;
Jac(2,:) = 2*A/W*GradA-A^2/W^2*GradW;
