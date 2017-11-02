A = importdata('A.txt');
[m, m] = size(A);
b = zeros(m, 1);

b(1) = 1;
A(1,:) = 1;

for i = 1:m
    A(i,:) = A(i,:) / A(i,i);
end

n = 1;

C11 = A(1:n, 1:n);
C12 = A(1:n, n+1:m);
A22 = A(n+1:m, n+1:m);
A21 = A(n+1:m, 1:n);

C21 = linsolve(A22, A21);

CD = C11 - C12 * C21;


x1 = linsolve(CD, b(1:n));
x2 = - C21 * x1;

x = [x1', x2'];
x = x';
