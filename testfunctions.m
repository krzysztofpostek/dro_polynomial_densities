function [falpha,Alphaf,n] = testfunctions(num)

% testfunctions

%%% INPUT: 
% num: number of test function f

%%% OUTPUT:
% n number of variables of f
% falpha coefficients of f in multivariate Chebyshev basis
% Alphaf matrix of multi-indices corresponding to falpha

if num == 1
    
%Booth 
n=2;
falpha = [250 250 800 -340 -380 574]';

Alphaf = [2 0 ; 0 2 ; 1 1 ; 1 0 ; 0 1 ; 0 0];

end

if num == 2
    
    %Matyas
    n=2;
    falpha = [13 13 -48 26]';

    Alphaf = [2 0 ; 0 2; 1 1; 0 0];

end


if num == 3

% Motzkin
% f_min = 0

n=2;

falpha = [4 4 4 4 20 16 16 13]';

Alphaf = [4 0 ; 4 2; 2 4; 0 4; 2 2; 2 0; 0 2; 0 0];

end


if num == 4
%Three hump camel

n=2;

falpha = [(5^6)/192 1625/4 58725/64 25 12.5 14525/24]';

Alphaf = [6 0 ; 4 0; 2 0; 1 1; 0 2; 0 0];

end

if num == 5

% Syblinsky-Tang n=2

n=2;

falpha = [(25^2)/16 (25^2)/16  225/4 225/4 25/2 25/2 (275*n)/16]' ;

Alphaf = [4 0  ; 0 4  ;  2 0 ; 0 2 ;  1 0  ; 0 1  ; 0 0];

end 

if num == 6

% Syblinsky-Tang n=3

n=3;

falpha = [(25^2)/16 (25^2)/16  (25^2)/16 225/4 225/4 225/4 25/2 25/2 25/2 (275*n)/16]' ;

Alphaf = [4 0 0 ; 0 4 0 ; 0 0 4 ; 2 0 0; 0 2 0; 0 0 2; 1 0 0 ; 0 1 0 ; 0 0 1; 0 0 0];

end

if num == 7
    
    %Rosenbrock n=2
    n=2;
    falpha = [  12.5*2.048^4  -100*2.048^3  (0.5 + 50*2.048^2)*2.048^2 50*2.048^2 -4.096  -100*2.048^3 (n-1)*(1 +2.048^2*(37.5*2.048^2+50.5))]';
    Alphaf = [4 0 ; 2 1 ; 2 0; 0 2; 1 0 ; 0 1; 0 0];
end    



if num == 8
    
    %Rosenbrock n=3
    n=3;
    falpha = [  12.5*2.048^4  -100*2.048^3  (0.5 + 50*2.048^2)*2.048^2 50*2.048^2 -4.096  -100*2.048^3 12.5*2.048^4  -100*2.048^3  (0.5 + 50*2.048^2)*2.048^2 50*2.048^2 -4.096  -100*2.048^3 (n-1)*(1 +2.048^2*(37.5*2.048^2+50.5))]';
    Alphaf = [4 0 0; 2 1 0; 2 0 0; 0 2 0; 1 0  0; 0 1 0; 0 4 0 ; 0 2 1 ; 0 2 0; 0 0 2; 0 1 0 ; 0 0 1; 0 0 0];
end 

end