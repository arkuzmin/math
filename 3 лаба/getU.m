function [U, mean_precision] = getU()
    a = 1;
    b = 1;
    N = 50;
    eps = 0.00001;
    maxIt = 1000;
    
    U = PuassonFunc(a, b, N, eps, maxIt);
    printPlot(a, b, N, U);
    
    % ��������� U c ������� ����� � 2 ����
    U_less_step = PuassonFunc(a, b, 2*N - 1, eps, maxIt);    

    % ��������� ��������
    precision = zeros(N, N);
    for i = 1:N
        for j = 1:N
            precision(i,j) = abs(U(i,j) - U_less_step(2*i-1, 2*j-1));
        end;
    end;
    mean_precision = mean(precision(:));

function [U] = PuassonFunc(a, b, N, eps, maxIt)
    hx = a / (N - 1);  
    hy = b / (N - 1);

    [y x] = meshgrid(0 : hy : b, 0 : hx : a);  

    U = zeros(N);    

    for i = 1:N
        U(1, i) = phi_0(y(1, i));    
        U(N, i) = phi_a(y(N, i));
        U(i, 1) = psi_0(x(i, 1));
        U(i, N) = psi_b(x(i, N));
    end    

    currentU = U; 

    % �������� ��������� w ���
    % ������ ������� ����������
    w = 2 / (1 + sin(pi/N));    
    
    x = 0 : hx : a;
    y = 0 : hy : b;
    for it = 1 : maxIt    
        for i = 2 : N - 1
            for j = 2 : N - 1
                u_2 = hy^2 / 2 / (hx^2 + hy^2) * (currentU(i-1, j) + U(i+1, j)) + ...
                      hx^2 / 2 / (hx^2 + hy^2) * (currentU(i, j-1) + U(i, j+1)) + ...
                      hx^2 * hy^2 / 2 / (hx^2 + hy^2) * f(x(i), y(j));
                currentU(i,j) = w * u_2 + (1 - w) * U(i,j);
            end;
        end;
        
        % ������� ��������
        % ������� n = norm(v, p) ��������� p-����� ������� v
        % || v ||p = sum(abs(v).^p)^1/p.
        difference_norm = norm(currentU - U, 1);
        U = currentU;
        if difference_norm < (2 - w) * eps   
            break;
        end;    
    end;

function res = phi_0(y)
   % res = sin(pi * y)^2;
   res = 2 * y*y;

function res = phi_a(y)
   % res = exp(sin(pi * y)) - 1;
   res = 2 * sin(pi * y);

function res = psi_0(x)
  %  res = x*(1-x);
 res = x - x*x;

function res = psi_b(x)    
  %  res = x * (1 - x) * exp(x);
  res = 2 - 2 * x;
  
    
function res = f(x, y)
    res = 5 + x + y;
   % res = (x - y)^2;
    
function printPlot(a, b, N, U)
    [Ygrid Xgrid] = meshgrid(0:b/(N-1):b, 0:a/(N-1):a);
    surf(Ygrid, Xgrid, U);   
    xlabel('Y', 'FontSize', 12);  
    ylabel('X', 'FontSize', 12);
