function [DT,ddt,S,dtij,A,res] = MCC(t,fs,seis,dt_trace,tpre,tpost,tf_weighted)

% Identify number of traces and samples
[ntrace,nsamp] = size(seis);

% Define maximum lag
maxlag = round((tpost - tpre)*fs);

% Multi-channel cross-correlation
N    = sum(1:(ntrace-1));
row  = zeros(N,2); % Coefficient matrix row-index
col  = zeros(N,2); % Coefficient matrix column-index
aij  = zeros(N,2); % Coefficient matrix value weight (+/- 1 unless weighted)
dtij = zeros(N,1); % Delay between the ith and jth stations
kk   = 0; % Counter
for ii = 1:(ntrace-1)
    % Arrival window for ith-trace
    wi = (t >= (dt_trace(ii) + tpre)) & (t <= (dt_trace(ii) + tpost));
    
    % The ith-trace
    si = seis(ii,:);
    
    for jj = (ii+1):ntrace
        % Update counter
        kk = kk + 1;
        
        % Arrival window for jth-trace
        wj = (t >= (dt_trace(jj) + tpre)) & (t <= (dt_trace(jj) + tpost));
        
        % The jth-trace
        sj = seis(jj,:);
        
        % Normalized cross-correlatation of windowed traces
        [r,lags] = xcorr(si.*wi,sj.*wj,maxlag,'coeff');
        % Identify maximum of cross-correlation
        wij  = max(r);
        keep = (r == wij);
        % Store delay time (Ti - Tj) and matrix coefficient
        if tf_weighted
            % Normalize by correlation coefficient
            dtij(kk)  = wij*mean(lags(keep))/fs;
            aij(kk,:) = wij*[1, -1];
        else
            % No data weighting
            dtij(kk)  = mean(lags(keep))/fs;
            aij(kk,:) = [1, -1];
        end
        
        % Indexing
        row(kk,1) = kk;
        row(kk,2) = kk;
        col(kk,1) = ii;
        col(kk,2) = jj;
        
    end
    
end
% Coefficient matrix
A = sparse(row(:),col(:),aij(:),N,ntrace);
% Add zero-mean constraint (do NOT use weighted mean)
A    = cat(1,A,ones(1,ntrace));
dtij = cat(1,dtij,0);

% Best-fit delay times
[DT,~,~,~]  = lsqr(A,dtij,[],1000);
res = dtij(:) - A*DT(:); % Residual vector for error estimation

% Estimate errors using Eq. 8 of Vandecar and Crosson, 1990.
% This is likely an underestimate!
ddt = zeros(ntrace,2);
for ii = 1:ntrace
    % 1) Error Estimate following Vandecar and Crosson (1990)
    rji = 0;
    rij = 0;
    % Sum of squared residuals associated with ith-trace
    for ji = 1:(ii-1)
        ind = sum((ntrace-(1:ii))) + (ji-ntrace);
        rji = rji + (res(ind)^2);
    end
    for ij = (ii+1):ntrace
        ind = sum((ntrace-(1:ii))) + (ij-ntrace);
        rij = rij + (res(ind)^2);
    end
    % Estimated standard error
    ddt(ii,1) = sqrt((1/(ntrace-2))*(rji + rij));
end


% Estimate errors based on similarity between traces and stacked waveform
% following Chevrot, GJI 2002.

% Define stacked trace
% Interpolate traces to align phases on current delays
seis = seis';
seis = interp2(1:ntrace,t(:),seis,repmat(1:ntrace,nsamp,1),...
    repmat(t(:),1,ntrace) + repmat(DT(:)',nsamp,1),'linear',0);
seis = seis';
% Stack traces. Always stack normalized traces to prevent biasing from any
% anomalous amplitudes.
w = (t >= tpre) & (t <= tpost);
% Amplitude normalized traces
seis = scale_traces(w,seis,1);
% Amplitude normalized stacked trace
S = mean(seis,1);

% Apply window to all seismograms
seis = repmat(w,ntrace,1).*seis;
S    = w.*S;

% Normalized auto-correlation of stack
AC = xcorr(S,S,'coeff');
% Compute correlation with stack
for ii = 1:ntrace
    % Normalized cross-correlation with stack
    [r,lags] = xcorr(S,seis(ii,:),'coeff');
    
    % Identify where the auto-correlation coefficient exceeds the
    % correlation coefficient for the ith-trace. The width of this region
    % defines the error (note auto-correlation is symmetric so we only need
    % to consider one intersection).
    r  = max(r);
    jn = find((AC >= r),1,'last');
    
    % Estimated error
    ddt(ii,2) = (((lags(jn+1)-lags(jn))/(AC(jn+1)-AC(jn)))*(r-AC(jn)) + lags(jn))/fs;
end
