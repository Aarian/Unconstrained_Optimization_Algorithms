function [] = problem1()
disp('--------------------------------------------------------------------------------------------------problem 1_a-------------------------------------------')
	cvx_begin
		variables x1 x2
		minimize(x1 +x2)
		subject to 
			-2*x1 -  x2 <= 0 
			-x1 - 3*x2 <= -1
			-x1 <= 0
			-x2 <= 0
		%norm (x,Inf) <= e
	cvx_end
	x1
	x2
	disp('##########################################################################################')	
	disp('----------------------------------------------------------------------------------------------problem 1_b-------------------------------------------')
	cvx_begin
		variables x1 x2
		minimize(norm ([ x1 ; x2 ],Inf))
		subject to 
			-2*x1 -  x2 <= 0 
			-x1 - 3*x2 <= -1
			-x1 <= 0
			-x2 <= 0
		%norm (x,Inf) <= e
	cvx_end
	x1
	x2
	disp('##########################################################################################')	
	disp('----------------------------------------------------------------------------------------------problem 1_c-------------------------------------------')
	cvx_begin
		variables x1 x2
		minimize(x1^2+9*x2^2)
		subject to 
			-2*x1 -  x2 <= 0 
			-x1 - 3*x2 <= -1
			-x1 <= 0
			-x2 <= 0
		%norm (x,Inf) <= e
	cvx_end
	x1
	x2
	
end