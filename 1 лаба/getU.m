function [U, N] = getU(l, T, p, t0, Q, c, a, b, sigma, u0, h, tau)
    format long;
    function W = W(t) 
        if (0 <= t && t < t0)
            W = 2*Q*(t0-t);
        else
        if (t >= t0)
            W = 0;
        end
        end
    end

    function K = K(u)    
        K = a + b * u^sigma;
    end
    
    M = l / h + 1;
    N = T / tau + 1;  
    
    x = 0 : h : l;
    t = 0 : tau : T;
    U = zeros(M, N);

    U(1, :) = u0;
    
    alpha = tau / (p*c*h^2);
    for j = 1 : M - 1
        U(j + 1, :) = U(j, :);
    
        A = zeros(N, N);
        F = zeros(N, 1);

        A(1, 1) = -1;
        A(1, 2) =  1;

        F(1) = 0;

        for i = 2 : N - 1
            Kminus = ( K(U(j+1, i - 1)) + K(U(j+1, i)) ) / 2;
            Kplus = ( K(U(j+1, i)) + K(U(j+1, i + 1)) ) / 2;

            A(i, i - 1) = -Kminus * alpha;
            A(i, i)     = ((Kminus + Kplus) * alpha + 1);
            A(i, i + 1) = -Kplus * alpha;

            F(i) = U(j, i);
        end

        A(N, N-1) = -1;
        A(N, N) = 1;
        Kn_ = ( K(U(j+1, N)) + K(U(j+1, N-1)) ) / 2;
        F(N) = W(t(j+1)) * h / Kn_;     
        U(j+1, :) = sweep_method(A, F); % A \ F;
    end
    
    % floor(X) возвращает значения, 
    % округленные до ближайшего целого<= X; 
    [XG TG] = meshgrid(t, x);
    surf(XG, TG, U); 
    ylabel('x', 'FontSize', 12);
    xlabel('t', 'FontSize', 12);  
end


