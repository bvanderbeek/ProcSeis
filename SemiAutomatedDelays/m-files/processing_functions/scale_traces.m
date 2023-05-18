function [S,A] = scale_traces(w,S,ichan)
% Scaling factor
if isempty(ichan)
    % Three-component amplitude scaling
    A = sqrt(sum(S.^2,3));
    A = mean(A(:,w),2);
else
    % Root-mean-squared channel scaling
    A = rms(S(:,w,ichan),2);
end
nsamp = size(S,2);
nchan = size(S,3);
S = S./repmat(A,1,nsamp,nchan);