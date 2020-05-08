%diary('D:\MSC\Boyd convex optimization\MZLGH\Codes\report\1\1a,b,c.txt')
disp('--------------------------------------------problem 1-------------------------------------------')
%problem1()









diary off
%diary on 
%diary('D:\MSC\Boyd convex optimization\MZLGH\Codes\report\2\2.txt')
disp('************************************************************************************')
disp('--------------------------------------------problem 2-------------------------------------------')

%problem2()




diary off
%diary on
%diary('D:\MSC\Boyd convex optimization\MZLGH\Codes\report\3\3.txt')
disp('************************************************************************************')
disp('--------------------------------------------problem 3-------------------------------------------')
%problem3()

diary off
disp('************************************************************************************')
disp('--------------------------------------------problem 4 SteepestDecent BackTrack-------------------------------------------')
%F1 = steepest_backtrack([-1,1]',.7,.7,500);
%F1 ;
%F2 = steepest_backtrack([2,1]',.7,.7,500);
%F2;
%F3 = steepest_backtrack([0,1]',.7,.7,500);
%F3;
disp('************************************************************************************')
disp('--------------------------------------------problem 4 Newton BackTrack -------------------------------------------')
%N1 = Newton_backtrack([-1,1]',.5,.1,500);
%N1; 
%N2 = Newton_backtrack([2,1]',.5,.1,500);
%N2;
%N3 = Newton_backtrack([0,1]',.7,.1,500);
%N3;

%diary('D:\MSC\Boyd convex optimization\MZLGH\Codes\report\4\SteepestDecent Cubic[0 1].txt')
disp('************************************************************************************')
disp('--------------------------------------------problem 4 SteepestDecent Cubic -------------------------------------------')

%steepest_cubic(x1,c,itr_num)
%disp ('initial point [-1 1]')
%F1 = steepest_cubic([-1,1]',.7,10);
%disp ('F*')
%F1(end) 
%disp ('initial point [2 1]')
%F2 = steepest_cubic([2,1]',.7,60);
%disp ('F*')
%F2(end)
%disp (' initial point [0 1]')
%F3 = steepest_cubic([0,1]',.1,10);
%disp ('F*')
%F3 (end)

disp('************************************************************************************')
disp('--------------------------------------------problem 4 Newton Cubic -------------------------------------------')
%disp ('initial point [-1 1]')
%F1 = NewtonCubic([-1,1]',.7,10);
%disp ('F*')
%F1(end) 
%disp ('initial point [2 1]')
%F2 = NewtonCubic([2,1]',.7,60);
%disp ('F*')
%F2(end)
%disp (' initial point [0 1]')
%F3 = NewtonCubic([0,1]',.1,10);
%disp ('F*')
%F3 (end)



diary off
diary('D:\MSC\Boyd convex optimization\MZLGH\Codes\report\5a\5a.txt')
disp('************************************************************************************')
disp('--------------------------------------------problem 5_A Trust-------------------------------------------')
disp('initial point = [2,2]') 
disp('Delta = 2')
X_star = TR([2,2]',1/32 , 5 ,2, 30)
disp('Delta = 1')
X_star = TR([2,2]',1/32 , 5 ,1, 30)
disp('Delta = 0.5')
X_star = TR([2,2]',1/32 , 5 ,.5, 30)

diary off
disp('************************************************************************************')
disp('--------------------------------------------problem 5_B DogLeg and Cauchy -------------------------------------------')
%disp('initial point = [2,2]')
%disp('\delta = 0.2')
%P_dog = dogleg ( [2,2]', .2) ; 
%P_dog
%disp('\delta = 1') 
%P_dog =dogleg ( [2,2]', 1) ;
%P_dog
%disp('\delta = 2')
%P_dog =dogleg ( [2,2]', 2) ;
%P_dog


%disp('\delta = 0.2')
%P_c = Cuachy ( [2,2]', .2);
%P_c
%disp('\delta = 1') 
%P_c = Cuachy ( [2,2]', 1) ;
%P_c
%disp('\delta = 2')
%P_c = Cuachy ( [2,2]', 2) ;
%P_c
%diary('D:\MSC\Boyd convex optimization\MZLGH\Codes\report\6\6.txt')

disp('************************************************************************************')
disp('--------------------------------------------problem 6 A -------------------------------------------')
%steepest_backtrack_six([1,1,1]',.7,.7,500)
%X_est_steepest = steepest_backtrack_six([1,1,1]', .7, .7, 500);
%X_est_steepest
disp('************************************************************************************')
disp('--------------------------------------------problem 6 B -------------------------------------------')
%ConjugateG([1 1 1]' , 3)


disp('************************************************************************************')
disp('--------------------------------------------problem 6 C -------------------------------------------')
%X_real = inv([3 -1 0 ; -1 3 -1 ; 0 -1 3] ) *  [1;2;3]


