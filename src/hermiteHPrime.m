function y = hermiteHPrime(n,x)
%HERMITEHPRIME - Returns the derivative of the nth Hermite polynomial
%   
%   Syntax
%       y = HERMITEHPRIME(n,x)

    arguments
        n array {mustBeInteger}
        x array {mustBeFinite}
    end

    if n == 0
        y = zeros(size(x));
    else
        y = 2*n .* hermiteH(n-1,x);
    end

end