function X_est = ConjugateG(x,iteration )
	A = [3 -1 0 ; -1 3 -1 ; 0 -1 3] ;
	b= [1;2;3] ; 
	r = ones(3,iteration) ; 
	p = ones(3,iteration) ;
	X = ones(3,iteration) ;
	alpha = ones(1,iteration) ;
	beta =  ones(1,iteration) ;
	r(: , 1) = A*x(:,1) -b ;
	p(:,1) = -r(: , 1) ;
	k = 1 ; 
	X(:,1) = x ; 
	while r(: , k) ~= 0
		alpha(k) = (- r(: , k)' * p(: , k))  / (p(:,k)' * A * p(:,k)) ;
		X (:,k+1) = X(:,k) + alpha(k) * p(:,k) ;
		r (:,k+1) = A*X(:,k+1) -b ;
		beta (k+1) = (r(:,k+1)' * A * p(:,k))/(p(:,k)' * A * p(:,k)) ; 
		p(:,k+1) = -r(:,k+1) + beta(k+1) * p(:,k) ;
		k = k+1 ;
	end
	X_est = X(:,end)
	X;	
	r
	
end

