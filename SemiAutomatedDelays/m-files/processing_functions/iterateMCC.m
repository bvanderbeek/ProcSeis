function [DT,ddt,n] = iterateMCC(t,S,fs,meanDelay,DT_i,tpre,tpost,tf_weighted,...
    tol,nmax)

% Iterative MCC
r = 1000*tol;
n = 0;
while (r > tol) && (n < nmax)
    % Call multi-channel cross-correlation
    [DT,ddt] = MCC(t - meanDelay,fs,S,DT_i,tpre,tpost,tf_weighted);
    
    % Update metrics
    r    = max(abs(DT - DT_i));
    DT_i = DT;
    
    % Update counter
    n = n + 1;
end

% Convergence check
if n >= nmax
    warning(['After ',num2str(n),' iterations, multi-channel cross-correlation failed to converge.']);
    fprintf(['\n Residual at convergence is: ',num2str(round(1000*r)),' ms. \n']);
else
    fprintf(['\n MCC converged in ',num2str(n),' iterations. \n']);
end
