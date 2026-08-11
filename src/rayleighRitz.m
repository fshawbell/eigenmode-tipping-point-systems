% Rayleigh–Ritz approximation of leading eigenvalue/eigenfunction
run setup.m

% Domain
a = -1; b = 1;    % adjust to your problem
nBasis = 5;      % number of trial basis functions

% Define p(x), q(x), w(x)
y = -0.5;
D = g(NaN).^2/2;
p = @(x) -D.*exp(U(x,y)/D);               
q = @(x) -U(x,y).*exp(U(x,y)/D);               
w = @(x) exp(U(x,y)/D);               

phi = @(k,x) sin(k*pi*(x - a)/(b - a));   % simple sine basis

% Build stiffness (A) and mass (B) matrices
A = zeros(nBasis);
B = zeros(nBasis);

for i = 1:nBasis
    for j = 1:nBasis
        dphi_i = @(x) (i*pi/(b - a))*cos(i*pi*(x - a)/(b - a));
        dphi_j = @(x) (j*pi/(b - a))*cos(j*pi*(x - a)/(b - a));
        fA = @(x) p(x).*dphi_i(x).*dphi_j(x) + q(x).*phi(i,x).*phi(j,x);

        fB = @(x) w(x).*phi(i,x).*phi(j,x);

        A(i,j) = integral(fA,a,b);
        B(i,j) = integral(fB,a,b);
    end
end

[vecs,vals] = eig(A,B);
lambda = diag(vals);
[lambdaSorted, idx] = sort(lambda,'descend');
phiCoeffs = vecs(:,idx(1));    % coefficients of leading eigenfunction

fprintf('Leading eigenvalue ≈ %g\n',lambdaSorted(1));

% Construct the approximate eigenfunction
xPlot = linspace(a,b,200);
yApprox = zeros(size(xPlot));
for k = 1:nBasis
    yApprox = yApprox + phiCoeffs(k)*phi(k,xPlot);
end

figure(1);
plot(xPlot,yApprox);
xlabel('x'); ylabel('Approx eigenfunction');
title(['Leading eigenfunction, \lambda ≈ ' num2str(lambdaSorted(1))]);


figure(2);
plot(xPlot,U(xPlot,y));
xlabel('x'); ylabel('U(x)');
