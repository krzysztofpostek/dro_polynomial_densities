function [bound] = Lasserre_bound_box(r,deg)

%%% INPUT: %%%
% 
% r: 2r is the max degree of the sos density

%%% OUTPUT: %%%
%  Lasserre upper bound of order r on f_min (K = [-1,1]) for f(x) = -x^deg

A = zeros(r+1,r+1) ;
B = A;

for i = 0:r
    for j = 0:r
        if even(i+j) == 1
            B(i+1,j+1) = 0.5/(i+j+1);
        end
        
        if even(i+j+deg) == 1
            A(i+1,j+1) = -0.5/(i+j+deg+1);
        end
    end
end    
            

       
       % solve generalised eig problem Ax = lambda Bx
       disp('solve generalised eig problem Ax = lambda Bx')
       
      [V,lambda] = eig(A,B);
      
      lambdamin = lambda(1);
      
bound = lambdamin;
    