function [A,lambda] = rayleighSchrodingerPert(G,l0,M)

    assert(ismatrix(G) && (size(G,1)==size(G,2)), 'G must be a square matrix.')
    assert(iscolumn(l0), 'l0 must be a column vector.')
    assert(size(l0,1)==size(G,2), 'Dimensions of l0 and G must match.')
    assert(isnumeric(M) && isscalar(M) && isfinite(M) && M == floor(M), 'M must be an integer scalar.')
    assert(M >= 1, 'M must be greater than or equal to 1.')

    N = size(l0,1);
    lambda = zeros(N,M+1);
    A = zeros(N,N,M+1);
    
    lambda(:,1) = l0;
    A(:,:,1) = eye(N);
    denom = l0.' - l0;

    for m = 2:M+1

        % Compute m-th corrections of eigenvalues
        lambda(:,m) = sum(A(:,:,m-1).*G.',2);
        
        % Compute m-th corrections of eigenfunctions
        if m == 2
            A(:,:,m) = (-A(:,:,m-1)*G) ./ denom;
        else
            extraTerms = sum(A(:,:,m-1:-1:2).*reshape(lambda(:,2:m-1),[N 1 m-2]),3);
            A(:,:,m) = (-A(:,:,m-1)*G  + extraTerms) ./ denom;
        end
        
        % Set diagonal to 0
        Am = A(:,:,m);     
        Am(eye(N) == 1) = 0; 
        A(:,:,m) = Am;
    
    end
end

