function F = steepest_cubic(xx1,c,itr_num) 

	f = @(x) 100*(x(2)-x(1)^2)^2+(1-x(1))^2 ;
	gf = @ (x) [-2*(1-x(1))-400*x(1)*(x(2)-x(1)^2) ; 200*(x(2)-x(1)^2)] ;
	phi = @ (Alpha,x,p) 100*((Alpha*p(2)+x(2))-(Alpha*p(1)+x(1))^2)^2+(1-(Alpha*p(1)+x(1)))^2;
	gphi = @ (Alpha,x,p) 200*(p(2) - 2*p(1)*(x(1) + Alpha*p(1)))*(x(2) + Alpha*p(2) - (x(1) + Alpha*p(1))^2) + 2*p(1)*(x(1) + Alpha*p(1) - 1);
	
	%d1 = @(alpha0,alpha1,x,p) gphi(alpha0) + gphi(alpha1) -3* (phi(alpha0,x,p) - phi(alpha1,x,p))/(alpha0-alpha1)
	%d2 =@ (alpha0,alpha1,x,p) sign(alpha1-alpha0) *sqrt(d1(alpha0,alpha1,x,p)^2 - gphi(alpha0)*gphi(alpha1))
	%d1 = @(alpha0,alpha1) subs(gphi,Alpha,alpha0) + subs(gphi,Alpha,alpha1) -3* (subs(phi,Alpha,alpha0) - subs(phi,Alpha,alpha1))/(alpha0-alpha1)
	%d2 =@ (alpha0,alpha1) sign(alpha1-alpha0) *sqrt(d1(alpha0,alpha1)^2 - subs(gphi,Alpha,alpha0)*subs(gphi,Alpha,alpha1))										
	%syms x2 x1 Alpha p1 p2 
	%%ee = p1+p2
	%phi= 100*((Alpha*p2+x2)-(Alpha*p1+x1)^2)^2+(1-(Alpha*p1+x1))^2;
	%gphi = diff(phi,Alpha)
	%x0 = [0,1]'
	alpha_arr  =  .9*ones(1,itr_num) ;
	X = ones(2,itr_num) ;
	F = ones(1,itr_num) ;
	%P=[0;0]
	x_nxt = xx1 ;
	%disp('GF')
	gf( [-1,1]' );
	for k = 1:itr_num
		X(:,k) = x_nxt ;
		P = -gf( X(:,k) ) ;
		ii=2;
		alpha_arr_check  =  .9*ones(1,10000);
		A_0 = 0.2 ;
		A_n1 = 4 ;
		A_p1 =1  ; 
		%disp('PHI')
		%disp(   eval( subs(phi, [Alpha , x1,x2,p1, p2 ] , [alpha_arr_check(ii) , X(1,k) ,X(2,k) , P(1) , P(2) ])) )
		%disp (  eval ( subs(phi, [Alpha , x1,x2,p1, p2 ] , [0 , X(1,k) ,X(2,k) , P(1) , P(2) ])) + c * alpha_arr_check(ii) * eval( subs(gphi, [Alpha , x1,x2,p1, p2 ] , [alpha_arr_check(ii) , X(1,k) ,X(2,k) , P(1) , P(2) ])) )
		while ( phi(A_0 , X(:,k) , P) > phi(0, X(:,k) , P) + c* A_0 *gphi(0 , X(:,k) , P)    )%eval( subs(phi, [Alpha , x1,x2,p1, p2 ] , [alpha_arr_check(ii) , X(1,k) ,X(2,k) , P(1) , P(2) ]))   <= eval( subs(phi, [Alpha , x1,x2,p1, p2 ] , [0 , X(1,k) ,X(2,k) , P(1) , P(2) ])) + c * alpha_arr_check(ii) * eval( subs(gphi, [Alpha , x1,x2,p1, p2 ] , [alpha_arr_check(ii) , X(1,k) ,X(2,k) , P(1) , P(2) ]))   )  
			
			%disp('in while')
			%disp(A_0)
			%disp(phi(A_0 , X(:,k) , P))
			%disp(phi(0, X(:,k) , P) + c* A_0 *gphi(0 , X(:,k) , P))
			
			A_0 = -.5*gphi(0 , X(:,k) , P)*A_n1 ^2/(phi(A_n1 , X(:,k) , P)-phi(0, X(:,k) , P)-gphi(0 , X(:,k) , P)*A_n1);
			
			%d1 = eval(subs(gphi,Alpha,alpha_arr_check(ii-1)) + subs(gphi,Alpha,alpha_arr_check(ii)) -3* (subs(phi,Alpha,alpha_arr_check(ii-1)) - subs(phi,Alpha,alpha_arr_check(ii)))/(alpha_arr_check(ii-1)-alpha_arr_check(ii)));
            %%subs (d1,)
			%d2 = sign(alpha_arr_check(ii)-alpha_arr_check(ii-1)) *sqrt(d1^2 - eval(subs(gphi,Alpha,alpha_arr_check(ii-1))*subs(gphi,Alpha,alpha_arr_check(ii))))	;
			%alpha_arr_check(ii+1) = alpha_arr_check(ii) - (alpha_arr_check(ii) - alpha_arr_check(ii-1)) *( (eval(subs(gphi,Alpha,alpha_arr_check(ii))) + d2 - d1 /(eval( subs(gphi,Alpha,alpha_arr_check(ii)) - subs(phi,Alpha,alpha_arr_check(ii-1))) +2*d2 ) )) ; 
			%d1
			%d2
			
			
			d1 = gphi(A_n1 , X(:,k) , P) + gphi(A_0 , X(:,k) , P) - 3* (phi(A_n1 , X(:,k) , P) - phi(A_0 , X(:,k) , P))/(A_0 - A_n1);
			d2 = sign (A_0 - A_n1)  *sqrt(d1^2 - gphi(A_n1 , X(:,k) , P) * gphi(A_0 , X(:,k) , P));
			A_p1 = A_0 - (A_0 - A_n1) *( gphi(A_0 , X(:,k) , P) + d2 - d1) / (A_0 - A_n1 +2*d2);
			%A_p1 = A_0

			
			%AB = [A_n1^2 ,-A_0^2 ; -A_n1^3,A_0^3] * [phi(A_0 , X(:,k) , P)-phi(0, X(:,k) , P)-gphi(0 , X(:,k) , P)*A_0; 
			%phi(A_n1 , X(:,k) , P)-phi(0, X(:,k) , P)-gphi(0 , X(:,k) , P)*A_n1 ] *1/(A_n1^2 * A_0^2*(A_0 - A_n1))
			%aa =  AB(1) ; 
			%bb = AB(2) ;
			%A_p1 = (-bb + sqrt(bb^2  - 3*aa*gphi(0 , X(:,k) , P)))/(3*aa)
			%
			%ii = ii + 1;
			A_n1 = A_0 ;
			A_0 = A_p1;
			%if abs(A_0 - A_n1) <= 0.001 || A_0 <= 0.01
			%	A_0 = .5*A_n1
			%end			
		end 
		alpha_arr(k) = A_0;
		F(k) = f(X(:,k)) ;
		x_nxt = X(:,k) + alpha_arr(k) * P ; 
		P;
		alpha_arr(k) ;
		
	end
	
%	cvx_begin
%		variable x(2)
%		minimize(100*(x(2)-x(1)^2)^2+(1-x(1))^2 )
%		subject to 
% 
%	cvx_end
%	x	
	disp('X')
		X(:,end)
	disp('Alpha')
		alpha_arr(end)
	%scatter(x,y,sz,c,'filled')
	%figure 
	%
	%subplot(1,2,1);
	%plot(1:itr_num , alpha_arr)
    %
	%
	%
	%title(strcat(' x_0 = ',mat2str(x1)))
	%xlabel('iteration')
	%ylabel('Alpha(BackTraking)')
	%
	%
	%subplot(1,2,2);
	%plot(1:itr_num , log(F))
	%%surf(X(1,:) , X(2,:) , F)
	%title(strcat('Steepest Decent ',  ' x_0 = ',mat2str(x1) , ' F* = ' ,mat2str(F(end)) , ' X* = ' ,mat2str(X(:,end))  ))
	%xlabel('iteration')
	%ylabel('log(F_k*)')
    
	
	return 

end