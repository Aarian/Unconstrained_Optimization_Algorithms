function P_c  = Cuachy(xx, Delta) 
	f = @(x) (x(2)-x(1)^2)^2+(1-x(1))^2 ;
	gf = @(x) [ 2*x(1) - 4*x(1)*(- x(1)^2 + x(2)) - 2 ; - 2*x(1)^2 + 2*x(2)];
	hf = @(x) [ 2+12*x(1)^2-4*x(1) , -4*x(1) ; -4*x(1) , 2 ];
	
	BB = hf(xx) ;
	gg = gf (xx) ;
	
	Toe = @(g,B) (g'*B*g <=0) * 1 + (g'*B*g >0) * min( 1 , norm(g,2)^3/(Delta*g'*B*g) ) ;
	
	P_c = (-Toe(gg,BB) * Delta/norm(gg,2)) * gg;
	min_val_Cauchy_model  = f(xx)+ gg'*P_c



end