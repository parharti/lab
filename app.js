const express = require('express');
const app = express();
const port = 3000;

// Route for returning LCM MATLAB code
app.get('/lcm', (req, res) => {
    const code = `
clc;
clear all;
format short;

% Matlab Code of Least Cost Method (LCM)
% Input Information
%% Input Phase
Cost = [6 4 1 5; 8 9 2 7; 4 3 6 4];
A = [14 16 5];
B = [6 10 15 4];

%% To check unbalanced/balanced Problem
if sum(A) == sum(B)
    fprintf('Given Transportation Problem is Balanced \\n');
else
    fprintf('Given Transportation Problem is Unbalanced \\n');
    if sum(A) < sum(B)
        Cost(end+1, :) = zeros(1, length(B));
        A(end+1) = sum(B) - sum(A);
    elseif sum(B) < sum(A)
        Cost(:, end+1) = zeros(length(A), 1);
        B(end+1) = sum(A) - sum(B);
    end
end

ICost = Cost;
X = zeros(size(Cost));
%% Initialize allocation
[m, n] = size(Cost);
% Finding No. of rows and columns
BFS = m + n - 1;  % Total No. of BFS

%% Finding the cell(with minimum cost) for the allocations
for i = 1:size(Cost, 1)
    for j = 1:size(Cost, 2)
        hh = min(Cost(:));  % Finding minimum cost value
        [Row_index, Col_index] = find(hh == Cost);  % Finding position of minimum cost cell
        x11 = min(A(Row_index), B(Col_index));
        [Value, index] = max(x11);  % Find maximum allocation
        ii = Row_index(index);  % Identify Row Position
        jj = Col_index(index);  % Identify Column Position
        y11 = min(A(ii), B(jj));  % Find the value
        X(ii, jj) = y11;
        A(ii) = A(ii) - y11;
        B(jj) = B(jj) - y11;
        Cost(ii, jj) = inf;
    end
end

%% Print the initial BFS
fprintf('Initial BFS =\\n');
IBFS = array2table(X);
disp(IBFS);

%% Check for Degenerate and Non Degenerate
TotalBFS = length(nonzeros(X));
if TotalBFS == BFS
    fprintf('Initial BFS is Non-Degenerate \\n');
else
    fprintf('Initial BFS is Degenerate \\n');
end

%% Compute the Initial Transportation cost
InitialCost = sum(sum(ICost .* X));
fprintf('Initial BFS Cost is = %d \\n', InitialCost);
    `;
    res.type('text/plain').send(code);
});

// Route for returning Steepest Descent MATLAB code
app.get('/steep', (req, res) => {
    const code = `
syms x1 x2
f1 = x1^2 - x1*x2 + x2^2;
fx = inline(f1);
fobj = @(x)fx(x(:, 1), x(:, 2));
grad = gradient(f1);
G = inline(grad);
gradx = @(x)G(x(:, 1), x(:, 2));
H1 = hessian(f1);
Hx = inline(H1);
xo = [1 0.5];
N = 2;
tol = 0.05;
iter = 0;
X = [];
while norm(gradx(xo)) > tol && iter < N
    X = [X; xo];
    S = -gradx(xo);
    H = Hx(xo);
    lam = S' * S ./ (S' * H * S);
    xnew = xo + lam * S';
    xo = xnew;
    iter = iter + 1;
end
disp(fobj(xo));
    `;
    res.type('text/plain').send(code);
});

// Route for returning Primal Simplex MATLAB code
// Route for returning Primal Simplex MATLAB code
app.get('/simplex', (req, res) => {
    const code = `
% Simplex Method
%max z=2x1+5X2
%x1+4x2<=24
%3x1+1x2<=21
%x1+x2<=9
clc
clear all
format short
Noofvariables=2;
C=[2 5];
a=[1 4; 3 1; 1 1]
b=[24; 21; 9]
s=eye(size(a,1))
A=[a s b]
cost=zeros(1,size(A,2))
cost(1:Noofvariables)=C
bv= Noofvariables+1:1:size(A,2)-1
zjcj=cost(bv)*A-cost
zcj=[zjcj; A]
simptable=array2table(zcj);
simptable.Properties.VariableNames(1:size(zcj,2))={'x_1','x_2','s_1','s_2','s_3','sol'}
RUN=true;
while RUN
if any(zjcj<0); %check for (most) negative value
    fprintf(' the current BFS is not optimal \n')
   zc=zjcj(1:end-1);
   [Enter_val, pvt_col]= min(zc) 
   if all(A(:,pvt_col)<=0)
    error('LPP is Unbounded all enteries are <=0 in column %d',pvt_col);
   else
       sol=A(:,end)
       column=A(:,pvt_col)
       for i=1:size(A,1)
         if column(i)>0
            ratio(i)= sol(i)./column(i)
         else
            ratio(i)=inf
         end
       end
       [leaving_val, pvt_row]=min(ratio)
   end
bv(pvt_row)=pvt_col
pvt_key=A(pvt_row, pvt_col)
A(pvt_row,:)=A(pvt_row,:)./pvt_key
for i=1:size(A,1)
    if i~=pvt_row
        A(i,:)=A(i,:)-A(i, pvt_col).*A(pvt_row,:)
    end
end
    zjcj=zjcj-zjcj(pvt_col).*A(pvt_row,:)
    zcj=[zjcj;A]
    table=array2table(zcj)
    table.Properties.VariableNames(1:size(zcj,2))={'x_1','x_2','s_1','s_2','s_3','sol'}
else
    RUN=false;
    fprintf('The current BFS is optimal \n')
end
end
`;
    res.type('text/plain').send(code);
});


// Route for returning Dual Simplex MATLAB code
// Route for returning Dual Simplex MATLAB code
app.get('/dual', (req, res) => {
    const code = `
format short
clc
clear all

C=[-2 0 -1 0 0 0]
A=[-1 -1 1 1 0 -5; -1 2 -4 0 1 -8]
ib=[4 5]
zjcj=C(ib)*A-C
RUN=true;
while RUN
    if any(A(:,size(A,2))<0)
        fprintf('the current BFS is not feasible')
        [lvg_val, pvt_row]=min(A(:,size(A,2)))
for i=1:size(A,2)-1
    if A(pvt_row,i)<0
        m(i)=zjcj(i)/A(pvt_row,i)
    else
         m(i)=-inf
     end
end
[ent_val, pvt_col]=max(m)
A(pvt_row,:)=A(pvt_row,:)/A(pvt_row,pvt_col)
for i=1:size(A,1)
     if i~=pvt_row
         A(i,:)=A(i,:)-A(i,pvt_col).*A(pvt_row,:)
     end 
end
ib(pvt_row)=pvt_col;
zjcj=zjcj-zjcj(pvt_col).*A(pvt_row,:)
ZCj=[zjcj;A]
        TABLE=array2table(ZCj);
        TABLE.Properties.VariableNames(1:size(ZCj,2))={'x_1','x_2','x3','s_1','s_2','sol'}
else
    RUN=false;
 fprintf('    current BFS is Feasible and Optimal   \n')
end
end
`;
    res.type('text/plain').send(code);
});


// Start the server
app.listen(port, () => {
    console.log(`Server is running on http://localhost:${port}`);
});
