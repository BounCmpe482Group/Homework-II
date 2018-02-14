function X = qr_householder(X)
    [m,n] = size(X);
    
    % Add one more row to X, so that v_k vectors can be stored
    X = [X; zeros(1,n)];
    
    for k = 1:n
        x = X(k:m, k);
        
        %We need sign of 0 to be 1 instead of 0, see sign_ext
        v_k = sign_ext(x(1))*norm(x)*eye(size(x)) + x;
        v_k = v_k / norm(v_k);
        
        % Construct upper triangular matrix
        X(k:m, k:n) = X(k:m, k:n) - 2*v_k*(v_k'*X(k:m, k:n));
        
        % Store v_k's
        X(k+1:m+1,k) = v_k;
    end
end


function r = sign_ext(x)
% SIGN_EXT  extend the sign function so that 0 has sign 1.
% 
% r = sign_ext(x)   sign of x

    if x == 0
        r = 1;
    else
        r = sign(x);
    end
end