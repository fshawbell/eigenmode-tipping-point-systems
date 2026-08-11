function [lambda, psi, p0, p1] = boundaryPert(Uhandle, D, x, opts)

    arguments
        Uhandle
        D
        x
        opts.reflectingBC (1,1) logical = false
        opts.fhandle (1,1) function_handle
    end
    
    if opts.reflectingBC && isempty(opts.fhandle)
        error('fhandle must be provided when reflectingBC is true.');
    end

    p0Inner = exp(Uhandle(x)/D);
    p0Outer = exp(-Uhandle(x)/D).*cumtrapz(x,p0Inner);
    c0 = 1/trapz(x,p0Outer);
    p0 = c0*p0Outer;
    
    p0Int = flip(cumtrapz(flip(x),flip(p0)));
    
    c1Inner = exp(Uhandle(x)/D).*p0Int;
    c1Outer = exp(-Uhandle(x)/D).*cumtrapz(x,c1Inner);
    c1 = -trapz(x,c1Outer)/trapz(x,p0Outer);
    
    p1Inner = exp(Uhandle(x)/D).*(p0Int+c1);
    p1 = exp(-Uhandle(x)/D)./D.*cumtrapz(x,p1Inner);

    if opts.reflectingBC
        p0Prime = c0 + opts.fhandle(x).*p0/D;
        p1Prime = (c1 + p0Int + opts.fhandle(x).*p1)/D;
        lambda = -p0Prime(end)/p1Prime(end);
    else 
        lambda = -p0(end)/p1(end);
    end
    
    psi = p0+lambda*p1;

end