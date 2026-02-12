function [U] = solve_U_NewtonRaphson_simp_stokes(psi, ak, c, delta, Ret, zo, flag, max_iterations, tolerance, initial_guess)
    U = initial_guess;
    
    for iter = 1:max_iterations
        % Calculate current function value
        Re = U * delta * Ret;
        c_fs = 0.0288 * Re^(-1/5) * (1 + 577*Re^(-6/5))^(2/3);
        c_f = 2 * ((1/2*c_fs)^3 + flag*((1/0.4)*(log(delta/zo) - psi))^(-6))^(1/3);
        L = (1/pi) * (1 - c/U)^2*( ak^2/4 + ak^4/4 + ak^6*(45/64));
        L = L + 0.5*c_f;
        
        % Function F(U) = U - 1/sqrt(L)
        F = U - 1/sqrt(L);
        
        % Calculate derivative F'(U)
        % Derivative of c_fs with respect to U
        dRe_dU = delta * Ret;
        dc_fs_dRe = 0.0288 * (-1/5) * Re^(-6/5) * (1 + 577*Re^(-6/5))^(2/3) + ...
                    0.0288 * Re^(-1/5) * (2/3) * (1 + 577*Re^(-6/5))^(-1/3) * 577 * (-6/5) * Re^(-11/5);
        dc_fs_dU = dc_fs_dRe * dRe_dU;
        
        % Derivative of c_f with respect to U
        term = (1/2*c_fs)^3 + flag*((1/0.4)*(log(delta/zo) - psi))^(-6);
        dc_f_dc_fs = 2 * (1/3) * term^(-2/3) * 3 * (1/2)^3 * c_fs^2;
        dc_f_dU = dc_f_dc_fs * dc_fs_dU;
        
        % Derivative of L with respect to U
        dL_dU = (1/pi) * 2 * (1 - c/U) * (c/U^2)*( ak^2/4 + ak^4/4 + ak^6*(45/64)) + 0.5 * dc_f_dU;
        
        % Derivative of F(U) with respect to U
        dF_dU = 1 + 0.5 * L^(-3/2) * dL_dU;
        
        % Newton-Raphson update with optional damping
        delta_U = -F/dF_dU;
        
        % Apply damping if needed (optional - can help with convergence)
        alpha = 1.0;  % No damping by default (set < 1 if oscillating)
        U_new = U + alpha * delta_U;
        
        % Check convergence
        error = abs(U_new - U) / abs(U);
        
        if error < tolerance
            fprintf('Newton-Raphson converged after %d iterations (error = %.2e)\n', iter, error);
            U = U_new;
            return;
        end
        
        U = U_new;  % Update for next iteration
    end
end