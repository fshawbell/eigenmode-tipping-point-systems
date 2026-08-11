function y = hermiteHPrime(n,x)

    if n == 0
        y = zeros(size(x));
    else
        y = 2*n .* hermiteH(n-1,x);
    end

end