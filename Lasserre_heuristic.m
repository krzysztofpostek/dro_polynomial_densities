function [bound] = Lasserre_heuristic(R,r)

%%% INPUT: %%%
% R: order of the upper bound
% r: 2r is the max degree of the sos density

%%% OUTPUT: %%%
%  upper bound of order R on f_min (as defined in paper of Etienne, Krzysztof and Daniel.) 


% calls subroutine nmultichoosek.m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5

% get test function f of degree d

d = 1;
n=3;
fAlpha = [ 1 1 1]';

Alphaf = eye(n);

% n = number of variables of f
% falpha vector with nonzero coeficients of f in the canonical multivariate
% basis
% Alphaf matrix of monomials corresponding to falpha






  
% Generate multi-indices for canonical multivariate basis

disp('Build table of monomials')
Alpha =  [zeros(1,n)];

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
      lamda_min = lambda(1)
      x = x/sqrt(x'*B*x);
      
      if order<R
      
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

        newmoments = zeros(RowsAlpha,1);
        Bb = zeros(RowsAlpha,nchoosek(n+2*r,2*r));

        for i = 1:RowsAlpha
        for j = 1:nchoosek(n+2*r,2*r)

        alpha = Alpha(i,:);
        beta = Alpha(j,:);
        gamma = alpha + beta;
        if sum(gamma) <= 2*r*R+d
            Bb(i,j) = moments(findrows(Alpha,gamma));
        end  
    end
end
moments = Bb*h;
clear newmoments
      end
end
bound = lambda_min;
    