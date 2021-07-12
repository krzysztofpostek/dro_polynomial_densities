%addpath('C:\Users\65193kpo\Downloads\cvx-w64\cvx');
%cvx_setup
addpath(genpath('C:\Users\65193kpo\SOSTOOLS.301'));
cvx_solver sedumi;

M = 2;
LB = [0.8; 0.7]; % Lower bounds on asset returns
UB = [1.2; 1.3]; % Upper bounds on asset returns 
allocations = [0.75; 0.25]; % Asset allocations
VaR_threshold = 0.9; % We want to calculate P(allocations'*(LB + UB .* z) <= VaR_threshold)

% Input data - exponents of the monomials whose moments we know
I_degrees = [0 0;
            1 0;
            0 1];
        
% Input data - known moments of the monomials on the support under the
% polynomial density
I_moments = [1; 0.5 * ones(M,1)];

polynomial_degree_range = [5];
wc_probs_table = zeros(length(polynomial_degree_range), 1);
integral_TV_table = zeros(length(polynomial_degree_range), 1);
time_table = zeros(length(polynomial_degree_range), 3);

for i_polynomial_degree = 1:length(polynomial_degree_range)
    
    poly_degree = polynomial_degree_range(i_polynomial_degree);
    [A, B] = wc_probability_polynomial2_MC(LB, UB, allocations, VaR_threshold, I_degrees, I_moments, poly_degree);
    wc_probs_table(i_polynomial_degree) = A;
    
    subplot(round(sqrt(length(polynomial_degree_range))), length(polynomial_degree_range)/ round(sqrt(length(polynomial_degree_range))), i_polynomial_degree);
    
    if(M == 2)
        plot_the_polynomial(M, poly_degree, B);
    end
    
end

% fprintf('%3.0f & %6.4f & %6.4f & %6.4f & %6.4f \\\\ \n', [polynomial_degree_range' wc_probs_table MC_probs_table integral_TV_table time_table(:,1)]')