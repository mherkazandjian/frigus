% A0 original matrix
% As scaled matrix by column, the diagonal elements have been recomputed, 
%    after adding the normalization line the first row is scaled by the
%    maximum elements
% Ar scaled by row
% A   scaled by column
%
clc; clear;
A0 = importdata('../../data/matlab_dev_matrices/A.txt'); 
A = A0;
[m, m] = size(A);
b = zeros(m, 1);
d = (diag(A));


As=A;
for i = 1:m    
    As(i,i) = -sum(sort(A([1:i-1,i+1:m],i)));
end

ds = (diag(As));
d = (diag(A));

Ar=A;
Ar(1,:)=1;
br = zeros(m, 1); br(1)=1;

for i = 1:m
    As(:,i) = As(:,i) / As(i,i);
    A(:,i) = A(:,i) / A(i,i);
    Ar(i,:)=Ar(i,:) /Ar(i,i);
end


dmaxI=max(ds);
b(1) = dmaxI;
As(1,:) = dmaxI./ds;



n = 4;

C11 = As(1:n, 1:n);
C12 = As(1:n, n+1:m);
A22 = As(n+1:m, n+1:m);
A21 = As(n+1:m, 1:n);

C21 = linsolve(A22, A21);


CD = C11 - C12 * C21;


x1 = linsolve(CD, b(1:n));
x2 = - C21 * x1;

y = [x1', x2']';
x = diag(1./ds)*y;

yf=As\b;
xf = diag(1./ds)*yf;


xrf=Ar\br;
disp('Condition number of A22, Ar and As')
disp([cond(A22) cond(Ar) cond(As)])
disp('First two components of x, xrf and xf')
disp([x(1:2)  xrf(1:2) xf(1:2)])
