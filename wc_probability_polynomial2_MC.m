% This is code for computing the Value-at-Risk of a portfolio under
% polynomial stuff

function [worst_case_probability, density_coefficients] = wc_probability_polynomial2(LB, UB, allocations, VaR_threshold, I_degrees, I_moments, poly_degree)

    % Definining the symbolic variable with polynomial density
    M = size(I_degrees, 2);
    z = sym('z',[M 1]);

    % Building up the matrix with monomial exponents of the density polynomial
    Density_mon_exponents = monomials_exponents(z, [0:poly_degree]);
    global Density_mon_exponents_2r
    Density_mon_exponents_2r = monomials_exponents(z, [0:2 * poly_degree]);
    clear z;
    
    full(Density_mon_exponents_2r)'
    full(Density_mon_exponents)'

    % Building up objective function coefficient vector

    w = (UB - LB) .* allocations;
    r = VaR_threshold - allocations' * LB;
    
    Objective_coeffs = zeros(size(Density_mon_exponents_2r, 1), 1);
    
    N_sim = 10^5;
    z_sim = rand(N_sim, M);
    
    % Computing the objective coefficient vector
    
    for i_sim = 1:N_sim
        if(z_sim(i_sim,:) * w <= r)
            Objective_coeffs = Objective_coeffs + prod(repmat( z_sim(i_sim,:), [size(Density_mon_exponents_2r, 1) 1]) .^ Density_mon_exponents_2r , 2) / N_sim;
        end
    end
    
    % Building up the matrix with linear constraints on the coefficients of the
    % density polynomial

    Density_coefficients_LHS = zeros(size(I_degrees, 1), size(Density_mon_exponents_2r, 1));
    
    for i_sim = 1:N_sim
        for i_moment_condition = 1:size(I_degrees,1)
            Density_coefficients_LHS(i_moment_condition, :) = Density_coefficients_LHS(i_moment_condition, :) + ...
                prod(repmat( z_sim(i_sim,:), [size(Density_mon_exponents_2r, 1) 1]) .^ (Density_mon_exponents_2r + repmat(I_degrees(i_moment_condition,:), [size(Density_mon_exponents_2r, 1) 1])), 2)' / N_sim;
        end
    end

    cvx_begin
        variable poly_coeffs(size(Density_mon_exponents_2r, 1), 1);
        variable Q(size(Density_mon_exponents, 1), size(Density_mon_exponents, 1)) symmetric;
        expression poly_coeffs_constraint(size(Density_mon_exponents_2r, 1),1);
        maximize(poly_coeffs' * Objective_coeffs);
        subject to

            Density_coefficients_LHS * poly_coeffs == I_moments;

            for i_row = 1:size(Density_mon_exponents,1)
                for i_col = 1:i_row-1

                    degree_in_main_poly = Density_mon_exponents(i_row,:) + Density_mon_exponents(i_col,:);
                    main_poly_coeff_index = find(all(repmat(degree_in_main_poly, [size(Density_mon_exponents_2r,1) 1]) == Density_mon_exponents_2r, 2));
                    poly_coeffs_constraint(main_poly_coeff_index) = poly_coeffs_constraint(main_poly_coeff_index) + 2 * Q(i_row, i_col);
                    [Density_mon_exponents(i_row,:); Density_mon_exponents(i_col,:)];
                    Density_mon_exponents_2r(main_poly_coeff_index, :);

                end

                degree_in_main_poly = Density_mon_exponents(i_row,:) + Density_mon_exponents(i_row,:);
                main_poly_coeff_index = find(all(repmat(degree_in_main_poly, [size(Density_mon_exponents_2r,1) 1]) == Density_mon_exponents_2r, 2));
                poly_coeffs_constraint(main_poly_coeff_index) = poly_coeffs_constraint(main_poly_coeff_index) + Q(i_row, i_row);
            end

            poly_coeffs_constraint == poly_coeffs;
            Q == semidefinite(size(Density_mon_exponents,1));
    cvx_end
    
    double(poly_coeffs' * Objective_coeffs)
    
    worst_case_probability = cvx_optval;
    density_coefficients = double(poly_coeffs);
    
%     N_sim = 10^5;
%     z_sim = 10 * randn(N_sim, M);
%     nnegative = 0;
%     
%     for i_sim = 1:N_sim 
%         point_c = z_sim(i_sim, :) + [0.5 0.5];
%         
%         val_poly = sum(prod(repmat(point_c, [size(Density_mon_exponents_2r, 1) 1]) .^ Density_mon_exponents_2r, 2) .* poly_coeffs);
%         
%         if(val_poly < - 10^5)
%            nnegative = nnegative + 1; 
%         end
%     end
%     
%     nnegative
end

