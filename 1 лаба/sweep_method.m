function X = sweep_method(M, V)
clc
% MX=V, M - трехдиагональная матрица
%M =[1 2 0 0 0
%    3 5 3 0 0
%    0 2 3 3 0 
%    0 0 1 1 3   
%    0 0 0 1 1];     
%V=[2 4 3 3 4]';
 
A = diag(M,-1); % - поддиагональ матрицы коэффициентов
B = diag(M);    % - главная диагональ матрици коэффициентов
C = diag(M,1);  % - наддиагональ матрицы коэффициентов
D = V;          % - вектор правой части системы
n = length(B);
alpha = zeros(1,n);
betta = zeros(1,n);

alpha(1) = C(1) / B(1);
for j=1:n-2
    alpha(j+1) = C(j+1) / (B(j+1) - A(j)*alpha(j));
end

betta(1) = D(1) / B(1);
for j=1:n-1
    betta(j+1) = (D(j+1) - A(j)*betta(j)) / (B(j+1) - A(j)*alpha(j));
end
X(n) = betta(n);

for j=n-1:-1:1
    X(j) = betta(j) - alpha(j)*X(j+1);
end
X = X';
end