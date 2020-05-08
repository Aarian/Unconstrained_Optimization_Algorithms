function F = steepest_backtrack_six(x1,ro,c,itr_num)
	A = [3 -1 0 ; -1 3 -1 ; 0 -1 3] ;
	b= [1;2;3] ; 
	
	f = @(x) (3*x(1) - x(2) -1)^2 + (-x(1) + 3*x(2) - x(3) -2)^2 +(-x(2)+3*x(3) -3)^2;
	gf = @ (x) [20*x(1) - 12*x(2) + 2*x(3) - 2 ;22*x(2) - 12*x(1) - 12*x(3) - 4 ; 2*x(1) - 12*x(2) + 20*x(3) - 14] ;
	
	%x0 = [0,1]'
	alpha  =  .9*ones(1,itr_num) ;
	X = ones(3,itr_num) ;
	F = ones(1,itr_num) ;
	%P=[0;0]
	x_nxt = x1 ;
	for k = 1:itr_num
		X(:,k) = x_nxt ; 
		P = -gf( X(:,k) ) ;
		while ( f( X(:,k) + alpha(k) * P) > f( X(:,k)) + c*alpha(k)* gf( X(:,k))' * P ) 
			alpha(k) = ro * alpha (k) ; 	
		end 
		F(k) = f(X(:,k)) ;
		x_nxt = X(:,k) + alpha(k) * P ; 
		
	end
	
%	cvx_begin
%		variable x(2)
%		minimize(100*(x(2)-x(1)^2)^2+(1-x(1))^2 )
%		subject to 
%
%	cvx_end
%	x	
	disp('X')
		F =X(:,end);



	
	return 

end