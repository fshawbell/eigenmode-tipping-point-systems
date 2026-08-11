function [t,x] = eulerMaruyama(x0,tMax,nx,d,fHandle,gHandle)
    Dt=tMax/nx; 
    x=zeros(d,nx+1);
    t=linspace(0,tMax,nx+1); 
    sqrtDt=sqrt(Dt);
    x(:,1)=x0; xi=x0;
    for i=1:nx 
        ti = t(i);
        dW=sqrtDt*randn(d,1);
        xiNext=xi+Dt*fHandle(xi,ti)+gHandle(xi,ti).*dW;
        x(:,i+1)=xiNext; 
        xi=xiNext;
    end
end