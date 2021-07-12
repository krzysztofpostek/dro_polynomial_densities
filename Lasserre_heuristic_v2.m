function [bound] = Lasserre_heuristic_v2(r,R)

%%% INPUT: %%%
% R: order of the upper bound
% r: 2r is the max degree of the sos density

%%% OUTPUT: %%%
%  upper bound of order R on f_min (as defined in paper of Etienne, Krzysztof and Daniel.) 


% calls subroutine nmultichoosek.m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5

% get test function f of degree d

% n = number of variables of f
% falpha vector with nonzero coeficients of f in the canonical multivariate
% basis
% Alphaf matrix of monomials corresponding to falpha


%Matyas function: f(x,y) = 26(x^2+y^2) -48xy
% d = 2;
% n=2;
% fAlpha = [ 26 26 -48]';
% 
% Alphaf = [2 0; 0 2 ; 1 1];



% Motzkin polynomial f(x,y) = 64(x^4y^2 + x^2y^4)-48x^2y^2, f_min = -1
% d = 6;
% n=2;
% fAlpha = [1 64 64 -48]';
% 
% Alphaf = [0 0 ; 4 2; 2 4 ; 2 2];

% %x^2-y^2
 d=2;
 n=2;
 fAlpha = [-1 -1]';
 Alphaf = [2 0 ; 0 2];




  
% Generate multi-indices for canonical multivariate basis

disp('Build table of monomials')
Alpha =  [zeros(1,n)];
%Alpha will be a table of monomials in the canonical basis, with each row giving the exponents of a monomial 
for m = 1:2*r*R+d

   combs=nmultichoosek(1:n, m);


   Alpham = zeros(length(combs),n);
   for i = 1:n
       BB = combs;
       BB(BB~=i) = 0;
       BB = BB/i;
       Alpham(:,i) = sum(BB,2);
   end
   Alpha = [Alpha ; Alpham] ;

    % construct the moments for the first measure
   
end

clear combs Alpham BB
RowsAlpha = length(Alpha);

% construct the moments for the initial (Lebesgue) measure on [-1,1]^n
disp('Constructing initial moment sequence')

moments = zeros(RowsAlpha,1);
for i = 1:RowsAlpha
    if even(Alpha(i,:)) == 1
      moments(i) = (2^n)./(prod((Alpha(i,:)+1)))   ;
    end
end    

for order = 1:R
   

% form the matrices A and B of generalised eig problem
disp('Form the matrices A and B of generalised eig problem')

% form matrix B 
B = zeros(nchoosek(n+r,r),nchoosek(n+r,r));

for i = 1:nchoosek(n+r,r)
    for j = 1:nchoosek(n+r,r)
        
        alpha = Alpha(i,:);
        beta = Alpha(j,:);
        gamma = alpha + beta;
        B(i,j) = moments(findrows(Alpha,gamma));

    end
end

A = zeros(nchoosek(n+r,r),nchoosek(n+r,r));
for k = 1:length(Alphaf)       
       for i = 1:nchoosek(n+r,r)
            for j = 1:nchoosek(n+r,r)
                alpha = Alpha(i,:);
                beta = Alpha(j,:);
                gamma = alpha + beta + Alphaf(k,:);
                A(i,j) = A(i,j) + fAlpha(k)*moments(findrows(Alpha,gamma));
            end
       end
end      
       
       % solve generalised eig problem Ax = lambda Bx
       disp('solve generalised eig problem Ax = lambda Bx')
       
      [V,lambda] = eig(A,B);
      x = V(:,1);
      lambdamin = lambda(1)
      x = x/sqrt(x'*B*x);
      
      
      
        % form optimal density
        disp('Forming optimal density h')
        h = zeros(nchoosek(n+2*r,2*r),1);
        for i = 1:nchoosek(n+r,r)
            for j = 1:nchoosek(n+r,r)
                alpha = Alpha(i,:);
                beta = Alpha(j,:);
                gamma = alpha + beta ;
                h(findrows(Alpha,gamma)) = h(findrows(Alpha,gamma)) + x(i)*x(j);
            end
        end



        % update moments
        disp('Update moments')

        %compute the number of moments that need to be updated
        rowlength = 1;
        for i =1:2*r*(R-order)+d
            rowlength= rowlength + length(nmultichoosek(1:n, i));
        end    
        
        
        Bb = zeros(rowlength,nchoosek(n+2*r,2*r));

        for i = 1:rowlength
            for j = 1:nchoosek(n+2*r,2*r)

            alpha = Alpha(i,:);
            beta = Alpha(j,:);
            gamma = alpha + beta;
            Bb(i,j) = moments(findrows(Alpha,gamma));

            end
        end
        moments = Bb*h;

      
end
bound = lambdamin;
    