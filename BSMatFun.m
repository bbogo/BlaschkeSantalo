function val = Fun(x,opt)
	switch opt.dim
	case 3
		tr = x(1)+x(2);
		dt = x(1)*x(2)-x(3)^2;
		val = [tr,dt];
	case 6
		tr = x(1)+x(2)+x(3);
		dt = x(1)*x(2)*x(3)+2*x(4)*x(5)*x(6)-x(1)*x(6)^2-x(2)*x(5)^2-x(3)*x(4)^2;
		val = [tr,dt];
	case 10
		tr = x(1)+x(2)+x(3)+x(4);
		A = [x(1) x(5) x(8) x(10);
		     x(5) x(2) x(6) x(9);
		     x(8) x(6) x(3) x(7);
		     x(10) x(9) x(7) x(4)];    
		dt = det(A);
		val = [tr,dt];
	end
