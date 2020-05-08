function X_star = TR (xx,eta, Delta_hat,Delta0,itration)
	f = @(x) (x(2)-x(1)^2)^2+(1-x(1))^2 ;
	gf = @(x) [ 2*x(1) - 4*x(1)*(- x(1)^2 + x(2)) - 2 ; - 2*x(1)^2 + 2*x(2)];
	hf = @(x) [ 2+12*x(1)^2-4*x(1) , -4*x(1) ; -4*x(1) , 2 ];
	g= gf(xx);
	B = hf(xx) ;
	mdl = @ (p) f(xx) + g'*p + p'* B *p;
	
	%gf(xx) ;
	%hf(xx) ;
	%P_dog = dogleg ( xx, 10)
	
	X = ones(2,itration) ;
	P_dog = ones(2,itration) ;
	Delta = ones(1,itration) ;
	Delta(1) = Delta0 ;
	X(:,1) = xx; 
	for k = 1:itration 
		P_dog(:,k) = dogleg ( X(:,k), Delta(k)) ;
		ro = (f(X(:,k)) - f(X(:,k)+P_dog(:,k)) ) / (mdl([0 0]') - mdl(P_dog(:,k)) );
		
		if ( ro <.25 )
			Delta(k+1) = Delta(k)/4 ;
		else
			if (ro>.25  & norm(P_dog(:,k),2)==Delta(k))
				Delta(k+1) = min(2*Delta(k) , Delta_hat) ; 
			else 
				Delta(k+1) = Delta(k) ;
			end
		end
		if(ro > eta) 
			X(:,k+1) = X(:,k) +P_dog(:,k) ;
		else
			X(:,k+1) = X(:,k) ;
		end
	
	
	end
	X_star = X(:,end);
	X ;
	
end