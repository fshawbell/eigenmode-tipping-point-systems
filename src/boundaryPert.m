function [lambda, psi, p0, p1] = boundaryPert(Uhandle, D, xmesh, opts)
%BOUNDARYPERT - Approximated leading eigenmode of the Fokker-Planck operator
%   This MATLAB function returns the approximated leading eigenvalue and
%   eigenfunction of the Fokker-Planck (FP) operator on a discretised 
%   domain. 
%   
%   Syntax
%       [lambda, psi] = BOUNDARYPERT(Uhandle, D, x)
%       [lambda, psi] = BOUNDARYPERT(Uhandle, D, x, opts)
%       [lambda, psi, p0, p1] = BOUNDARYPERT(Uhandle, D, x, opts)
%
%   Input Arguments
%       Uhandle - FP potential function
%           function_handle
%       D - FP diffusion constant
%           scalar
%       xmesh - Spatial mesh
%           vector
%       UPrimehandle - Derivative of FP potential function
%           function_handle
%       neumannBC = false - Right hand boundary condition
%           logical
%   
%   Output Arguments
%       lambda - Approximated eigenvalue
%           scalar
%       psi - Approximated eigenfunction
%           vector
%       p0 - Unperturbed solution
%           vector
%       p1 - First order correction
%           vector

    arguments
        Uhandle function_handle
        D (1,1) {mustBeScalarOrEmpty,mustBeNonempty}
        xmesh {mustBeVector,mustBeNumeric,mustBeReal}
        opts.UPrimehandle (1,1) function_handle
        opts.neumannBC (1,1) logical = false
    end
    
    if opts.neumannBC && isempty(opts.UPrimehandle)
        error('UPrimehandle must be provided when neumannBC is true.');
    end

    p0Inner = exp(Uhandle(xmesh)/D);
    p0Outer = exp(-Uhandle(xmesh)/D).*cumtrapz(xmesh,p0Inner);
    c0 = 1/trapz(xmesh,p0Outer);
    p0 = c0*p0Outer;
    
    p0Int = flip(cumtrapz(flip(xmesh),flip(p0)));
    
    c1Inner = exp(Uhandle(xmesh)/D).*p0Int;
    c1Outer = exp(-Uhandle(xmesh)/D).*cumtrapz(xmesh,c1Inner);
    c1 = -trapz(xmesh,c1Outer)/trapz(xmesh,p0Outer);
    
    p1Inner = exp(Uhandle(xmesh)/D).*(p0Int+c1);
    p1 = exp(-Uhandle(xmesh)/D)./D.*cumtrapz(xmesh,p1Inner);

    if opts.neumannBC
        p0Prime = c0 + opts.UPrimehandle(xmesh).*p0/D;
        p1Prime = (c1 + p0Int + opts.UPrimehandle(xmesh).*p1)/D;
        lambda = -p0Prime(end)/p1Prime(end);
    else 
        lambda = -p0(end)/p1(end);
    end
    
    psi = p0+lambda*p1;

end