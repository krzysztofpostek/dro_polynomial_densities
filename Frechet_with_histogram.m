% Computing the worst-case probability that sum(z) >= 15, subject to
% knowledge on the probability of each z being in bins [0, 0.25], ...
% [4.75, 5] - these probabilities are set to be the ones from the marginal
% lognormal densities

% Uses function 
% adsimp() of Alan Genz
% cvx with some solver

% Set boxes for the distributions
bins = repmat([0:0.25:5]', [1 3]); %Setting the bins endpoints
m = [-0.3 0.4 0.8]; % Lognormal density means
v = [0.8 0.5 0.5]; % Lognormal density standard deviaation
b = 15; % Threshold
probabilities = zeros(size(bins, 1) - 1, 3); 
poly_degree = 4; % Polynomial degree
vertices = b * eye(3, 4); % Vertices of the simplex we integrate over

% Computing the real probabilities for the bins from the lognormal
% distribution function
for i_distribution = 1:3
    for i_bin = 1:size(bins, 1) - 1
        probabilities(i_bin, i_distribution) = logncdf(bins(i_bin + 1, i_distribution), m(i_distribution), v(i_distribution)) ...
            - logncdf(bins(i_bin, i_distribution), m(i_distribution), v(i_distribution));
    end
end

% Integration of a monomial z_1^monomial_1 * z_2^monomial_2 *
% z_2^monomial_2 times lognormal densities over a set [a b] x [0, +Infty] 
% x [0, + Infty] is by Fubini theorem decomposable to a product of (i)
% integral of z_1^monomial_1 times lognormal density of z_1 over the set
% [a,b], (ii) integrals of z_2^monomial_2, z_2^monomial_2 times their
% respective lognormal densities over [0, +Infty] computing moments of the marginal
% densities

% In Computed_moments(i_bin, i_distribution, i_degree + 1) we store the
% integral of z_{i_distribution}^{i_degree} over the bin 'i_bin' (this can
% be computed exactly)
Computed_moments = zeros([size(probabilities) 11]);

for i_distribution = 1:3
    for i_bin = 1:size(probabilities, 1)
        for i_degree = 0:10
            if(i_degree == 0)
                Computed_moments(i_bin, i_distribution, i_degree + 1) = probabilities(i_bin, i_distribution);
            else
                shifted_m = m(i_distribution) + i_degree * v(i_distribution)^2;
                
                if(i_bin == 1)
                    Computed_moments(i_bin, i_distribution, i_degree + 1) = exp(m(i_distribution) * i_degree + (i_degree * v(i_distribution))^2 / 2) * ...
                        normcdf(log(bins(i_bin + 1, i_distribution)), shifted_m, v(i_distribution));
                else
                    Computed_moments(i_bin, i_distribution, i_degree + 1) = exp(m(i_distribution) * i_degree + (i_degree * v(i_distribution))^2 / 2) * ...
                        (normcdf(log(bins(i_bin + 1, i_distribution)), shifted_m, v(i_distribution)) - ...
                        normcdf(log(bins(i_bin, i_distribution)), shifted_m, v(i_distribution)));
                end
            end
        end
    end
end

% Building up the matrix with monomial exponents of the density polynomial
% density will be encoded by nxn matrix Q where n = nb of rows in
% Density_mon_exponents
Density_mon_exponents = zeros(1, 3);

for i_m = 1:poly_degree
    
    combs=nmultichoosek(1:3, i_m);
    
    Alpham = zeros(length(combs), 3);
    for i = 1:3
        BB = combs;
        BB(BB~=i) = 0;
        BB = BB/i;
        Alpham(:,i) = sum(BB,2);
    end
    
    Density_mon_exponents = [Density_mon_exponents ; Alpham] ;
    
end

% Preallocating matrices
Objective_coeffs_matrix = zeros(size(Density_mon_exponents, 1), size(Density_mon_exponents, 1)); % Matrix for the objective function - contains integrals of monomials times lognormal densities' product over a simplex
Density_coefficients_support = zeros(size(Density_mon_exponents, 1), size(Density_mon_exponents, 1)); % Matrix for the condition that the density integrates to 1 - contains integrals of monomials times lognormal densities' product over nonnegative orthant

% Filling Objective_coeffs_matrix
for i_row = 1:size(Density_mon_exponents, 1)
    for i_col = 1:i_row-1 % Below the diagonal
        [poly_degree i_row i_col]
        monomial = Density_mon_exponents(i_row, :) + Density_mon_exponents(i_col, :); % Monomial corresponding to the (i_row, i_col)-th entry in the matrix defining the density
        % Here using Genz's function adsimp(number of variables, vertices
        % of the simplex, function, some integration parameter (number of
        % function evaluations I guess), integration precision)
        
        % r - integral value, e - error estimate
        [r, e] = adsimp(3, vertices, 1, @(x)(x(1)^monomial(1) * x(2)^monomial(2) * x(3)^monomial(3) / (x(1) * x(2) * x(3) * prod(v) * (2 * pi)^(3/2)) ...
            * exp(- (log(x(1)) - m(1))^2 ./ (2 * v(1)^2)) * exp(- (log(x(2)) - m(2))^2 ./ (2 * v(2)^2)) * exp(- (log(x(3)) - m(3))^2 ./ (2 * v(3)^2))), 10^4, 1e-5); 
        % Integration of a monomial w.r.t. product of lognormal densities over a simplex
        Objective_coeffs_matrix(i_row, i_col) = r;
        Objective_coeffs_matrix(i_col, i_row) = Objective_coeffs_matrix(i_row, i_col);
    end
    
    %... and on the diagonal
    monomial = 2 * Density_mon_exponents(i_row, :);
    [r, e] = adsimp(3, vertices, 1, @(x)(x(1)^monomial(1) * x(2)^monomial(2) * x(3)^monomial(3) / (x(1) * x(2) * x(3) * prod(v) * (2 * pi)^(3/2)) ...
        * exp(- (log(x(1)) - m(1))^2 ./ (2 * v(1)^2)) * exp(- (log(x(2)) - m(2))^2 ./ (2 * v(2)^2)) * exp(- (log(x(3)) - m(3))^2 ./ (2 * v(3)^2))), 10^4, 1e-5);
    % Integration of a monomial w.r.t. product of lognormal densities over a simplex
    Objective_coeffs_matrix(i_row, i_row) = r;
end

% Filling Density_coefficients_support - this is easy, just computing the
% moments of the distribution over the entire nonnegative orthant
for i_row = 1:size(Density_mon_exponents, 1)
    for i_col = 1:i_row-1
        monomial = Density_mon_exponents(i_row, :) + Density_mon_exponents(i_col, :);
        Density_coefficients_support(i_row, i_col) = prod(exp(monomial .* m + 0.5 * monomial.^2 .* v.^2));
        Density_coefficients_support(i_col, i_row) = Density_coefficients_support(i_row, i_col);
    end
    
    monomial = 2 * Density_mon_exponents(i_row, :);
    Density_coefficients_support(i_row, i_row) = prod(exp(monomial .* m + 0.5 * monomial.^2 .* v.^2));
end

% Moments_coefficients_support is a 3D array that we will use to ensure
% that per bin the probability is the right one. The 3-rd dimension is
% number_of_bins x number of marginal distributions
Moments_coefficients_support = zeros(size(Density_mon_exponents, 1), size(Density_mon_exponents, 1), size(probabilities, 1) * 3);

% Computing the integrals of the monomials of the density over each bin of
% each marginal distribution
for i_bin = 1:size(probabilities, 1)
    for i_distribution = 1:3
        i_moment = (i_distribution - 1) * size(probabilities, 1) + i_bin; % Indexing all these bins
        
        for i_row = 1:size(Density_mon_exponents, 1)
            for i_col = 1:i_row-1
                monomial = Density_mon_exponents(i_row, :) + Density_mon_exponents(i_col, :);
                integral = prod(exp(monomial(~([1:3] == i_distribution)) .* m(~([1:3] == i_distribution)) + ... 
                    0.5 * monomial(~([1:3] == i_distribution)).^2 .* v(~([1:3] == i_distribution)).^2));
                    % Product of moments of the two distributions that are not 'binned' here
                integral = integral * Computed_moments(i_bin, i_distribution, monomial(i_distribution) + 1); 
                % Times product the integral of z^alpha * lognormal density for the binned distribution, taken from the pre-computed Computed_moments
                
                Moments_coefficients_support(i_row, i_col, i_moment) = integral; % Insert into the matrix
                Moments_coefficients_support(i_col, i_row, i_moment) = Moments_coefficients_support(i_row, i_col, i_moment); % Matrix is symmetric...
            end

            monomial = 2 * Density_mon_exponents(i_row, :);
            integral = prod(exp(monomial(~([1:3] == i_distribution)) .* m(~([1:3] == i_distribution)) +...
                0.5 * monomial(~([1:3] == i_distribution)).^2 .* v(~([1:3] == i_distribution)).^2)); 
                % Product of moments of the two distributions that are not 'binned' here
            integral = integral * Computed_moments(i_bin, i_distribution, monomial(i_distribution) + 1); 
            % Times product the integral of z^alpha * lognormal density for the binned distribution, taken from the pre-computed Computed_moments
            Moments_coefficients_support(i_row, i_row, i_moment) = integral; % Insert into the matrix
        end
    end
end

% Optimization problem
cvx_begin
    variable Q(size(Density_mon_exponents, 1), size(Density_mon_exponents, 1)) symmetric; % Matrix defining the density
    minimize(sum(sum(Q .* Objective_coeffs_matrix))) % We minimize the probability of being inside the simplex
    subject to
        sum(sum(Density_coefficients_support .* Q)) == 1; % Density sums up to 1
        Q == semidefinite(size(Density_mon_exponents, 1)); % PSD condition

        for i_bin = 1:size(probabilities, 1)
            for i_distribution = 1:3
                i_moment = (i_distribution - 1) * size(probabilities, 1) + i_bin; % Indexing all these bins
                trace(Q' * Moments_coefficients_support(:, :, i_moment)) == probabilities(i_bin, i_distribution); % Matching the per bin marginal probability
            end
        end
cvx_end