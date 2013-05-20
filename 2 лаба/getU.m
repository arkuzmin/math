function U = getU()        
    l = 1;
    T = 1;
    a = 1;  
    
    h = 0.02;
    tau = 0.02;

    M = l / h + 1;
    N = T / tau + 1;
    
    
    
    U = Ufunc(a, h, tau, M, N);
    
    x = 0:h:l;
    t = 0:tau:T;
    
    [Xgrid Tgrid] = meshgrid(t, x);
    surf(Xgrid, Tgrid, U);   

    xlabel('t','FontSize',12);  
    ylabel('x','FontSize',12);

function U = Ufunc(a, h, tau, M, N)
    U = zeros(M, N);   
    
    % Определение значения сеточной функции на 1 слое
    % по начальному положению струны
    curr_x = 0;
    for i = 1 : M
        U(i, 1) = phi(curr_x);
        curr_x = curr_x + h;
    end;
    
    % Определение значения сеточной функции на 2 слое
    % по начальным скоростям
    curr_x = h;
    for i = 2 : M-1
        U(i, 2) = tau * psi(curr_x) + U(i, 1);
        curr_x = curr_x + h;
    end;
    
    % Определение значения функции в граничных узлах
    % из граничных условий
    curr_t = tau;
    for j = 2 : N
        U(1, j) = mu(curr_t);
        U(M, j) = v(curr_t);
        curr_t = curr_t + tau;
    end;
    
    % Определение значения сеточной функции
    % во внутренних узлах j+1 слоя
    % по разностному уравнению с исп-ем двух предыдущих слоев
    alpha = a^2 * tau^2 / h^2;
    for j = 2 : N-1
        for i = 2 : M-1
            U(i,j+1) = 2*U(i, j) - U(i, j-1) + alpha * (U(i-1, j) - 2*U(i, j) + U(i+1, j));
        end;
    end;

function y = phi(x)
 %   y = (x + 1) * sin(pi*x / 2);
  y = (x + 0.5) * (x + 1);

function y = psi(x)
  % y = x*x;
 y = cos(x + 0.5);
   
function y = mu(t)
 %   y = 0.5*t;
 y = 0.5;

function y = v(t)
 %   y = 2 - t^2;
 y = 3 - 2*t;