function F = steepest_backtrack(x1,ro,c,itr_num)
	f = @(x) 100*(x(2)-x(1)^2)^2+(1-x(1))^2 ;
	gf = @ (x) [-2*(1-x(1))-400*x(1)*(x(2)-x(1)^2) ; 200*(x(2)-x(1)^2)] ;
	
	%x0 = [0,1]'
	alpha  =  .9*ones(1,itr_num) ;
	X = ones(2,itr_num) ;
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
		X(:,end)
	disp('alpha')
		alpha(end)
	%scatter(x,y,sz,c,'filled')
	figure 
	
	subplot(1,2,1);
	plot(1:itr_num , alpha)
	%indexmax = find(max(F) == F);
	%Fmax = F(indexmax);
	%strmin = ['F* = ',num2str(F(end))];
	%text(itr_num/2,Fmax,strmin,'HorizontalAlignment','left');
    %
	%strmin = ['X* = ',num2str(F(:,end))];
	%text(itr_num/2,Fmax/2,strmin,'HorizontalAlignment','left');
	
	
	title(strcat(' x_0 = ',mat2str(x1)))
	xlabel('iteration')
	ylabel('alpha(BackTraking)')
	
	
	subplot(1,2,2);
	plot(1:itr_num , log(F))
	%surf(X(1,:) , X(2,:) , F)
	title(strcat('Steepest Decent ',  ' x_0 = ',mat2str(x1) , ' F* = ' ,mat2str(F(end)) , ' X* = ' ,mat2str(X(:,end))  ))
	xlabel('iteration')
	ylabel('log(F_k*)')

	
	return 

end