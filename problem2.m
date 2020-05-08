function []= problem2()
	n=3;
	P = [13 12 -2 ; 12 17 6 ; -2 6 12];
	q = [-22 ; -14.5 ; 13];
	r = 1;
	cvx_begin

	variables x(n) 
	minimize .5*transpose(x) * P * x+transpose(q) * x + r
	subject to
	%x(1) <= 1
	%x(2) <= 1
	%x(3) <= 1
	%-x(1) <= 1
	%-x(2) <= 1
	%-x(3) <= 1
		-1<= x <= 1
	cvx_end
	x


end