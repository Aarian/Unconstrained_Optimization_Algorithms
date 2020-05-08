function P_star  = dogleg(xx, Delta) 
	f = @(x) (x(2)-x(1)^2)^2+(1-x(1))^2 ;
	gf = @(x) [ 2*x(1) - 4*x(1)*(- x(1)^2 + x(2)) - 2 ; - 2*x(1)^2 + 2*x(2)];
	hf = @(x) [ 2+12*x(1)^2-4*x(1) , -4*x(1) ; -4*x(1) , 2 ];
	
	B = hf(xx) ;
	g = gf (xx) ;


	P_B = -inv(B) *g;
	P_U = - ((g'*g) /(g'*B*g) ) * g;
	
	f = @(x) (x(2)-x(1)^2)^2+(1-x(1))^2 ;
	P_hat = @(t) ((0<=t)&(t<=1)) * t' *P_U +  ((1<t)&(t<=2)) * (P_U + (t'-1) *(P_B - P_U));
	mdl = @ (p) f(xx) + g'*p + p'* B *p;
	syms t  x
	%eqn = norm(P_U + (t-1)*(P_B-P_U),2) == Delta
	eqn = (P_U(1) + (t-1)*(P_B(1)-P_U(1))) ^2 + (P_U(2) + (t-1)*(P_B(2)-P_U(2))) ^2 == Delta^2;
	tt = eval(solve(eqn,t));
	
	
	%P_star = 0 ;
	
	if(norm(P_B,2) <=  Delta || tt(2) >=2) 
		%disp('xx+P_B')
		xx+P_B;
		v_star = f(xx+P_B) ;
		P_star = P_B ;
	
	else
		if (imag(tt(2)) ~= 0)
			tt(2) = 0 ;
		end
			if (  0 <=tt(2) & tt(2) < 1  )
				%disp ('0< < 1')
				P_star = tt(2) * P_U ;
				xx+P_star ;
			else 
				if (1 <=tt(2) & tt(2) < 2)
					P_star =  (P_U + (tt(2)-1) *(P_B - P_U));
					%disp ('1< < 2')
					xx+P_star ;
				end
			end
		%end
	end		
	%disp('out of trust')
	%T = 0.1:.1:2;
	%sz = size(T)
	%lst = [] ;
	%my_p = 0;
	%toe = 0.1;
	%min_val_dog_model  = 1000 ;
	%for kk = 0:.1:2
	%	p_t =  P_hat( kk )
	%	%kk
	%	%disp('norm p_t ')
	%	%norm(p_t,2) 
	%	if(norm(p_t,2) <=  Delta)
	%		disp('||p_t|| <=  Delta')
	%		if(min_val_dog_model  > mdl(p_t) )
	%			min_val_dog_model  = mdl(p_t);
	%			my_p = p_t;
	%			toe = kk ;
	%		end
	%		%mdl_out = mdl(p_t) ;
	%		%lst = [lst,mdl_out];
	%	end
	%end
	%[m,i] = min ( mdl_out ) ;
	%P_star = P_hat( T(i) ) ;
	%P_star = my_p ;
	
	
	%disp ('min_val_dog_model  : ')
	%T(i) 
	%m
	%min_val_dog_model  
	%toe	
	%end

end