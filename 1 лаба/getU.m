

function [U, N] = getU()
    
    c = 0.5;
    a = 1;
    b = 2;
    sigma = 0.25;
    u0 = 0.1;
    l = 1;
    T = 1;
    p = 1;
    h = 0.02;
    tau = 0.02;
    t0 = 0.5;
    Q = 10;

    format long;
    function W = W(t) 
        if (0 <= t && t < t0)
            W = 2*Q*(t0-t);
          % W = 2*Q*t;
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

       % F(1) = 0;
       
       K0plus = ( K(U(j+1,2)) + K(U(j+1, 1)) ) / 2;
       F(1) = -W(t(j+1)) * h / K0plus;  
       

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
      %  F(N) = W(t(j+1)) * h / Kn_;     
        F(N) = 0;
        U(j+1, :) = sweep_method(A, F); % A \ F;
    end
    
    % floor(X) возвращает значения, 
    % округленные до ближайшего целого<= X; 
    [XG TG] = meshgrid(t, x);
    surf(XG, TG, U); 
    ylabel('x', 'FontSize', 12);
    xlabel('t', 'FontSize', 12);  
end


