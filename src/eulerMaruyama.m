function [t,x] = eulerMaruyama(tspan,x0,d,fhandle,ghandle)
%EULERMARUYAMA - Approximate realisations of SDEs
%   This MATLAB function approximates realisations of random processes as
%   described by an SDE with drift term fhandle and diffusion term ghandle.
%   
%   Syntax
%       [t,x] = EULERMARUYAMA(tspan,x0,d,fhandle,ghandle)
%
%   Input Arguments
%       tspan - Time steps
%           vector
%       x0 - Initial condition(s) 
%           vector
%       d - Number of realisations
%           integer
%       fhandle - Drift term
%           function_handle
%       ghandle - Diffusion term
%           function_handle
%   
%   Output Arguments
%       t - Evaluation points
%           vector
%       x - Realisations
%           array
    
    arguments
        tspan {mustBeVector, mustBeNumeric, mustBeReal}
        x0 {mustBeVector, mustBeNumeric, mustBeReal}
        d {mustBeScalarOrEmpty, mustBeInteger}
        fhandle function_handle
        ghandle function_handle
    end

    t=tspan; 
    nt = length(tspan);
    x=zeros(d,nt+1);
    x(:,1)=x0; xi=x0;
    for i=1:nt 
        ti = t(i);
        dt = t(i+1)-t(i);
        dW=sqrt(dt)*randn(d,1);
        xiNext=xi+dt*fhandle(xi,ti)+ghandle(xi,ti).*dW;
        x(:,i+1)=xiNext; 
        xi=xiNext;
    end
end