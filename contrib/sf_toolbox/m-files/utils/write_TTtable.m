function TTtable = write_TTtable(theFile,aModel,rminmax,zminmax,N,M)
% WRITE_TTTABLE: Creates a travel-time table structure for a 1D model using
% the TauP toolbox. Use the table to interpolate between 1D travel-time
% predictions instead of making many TauP calls from Matlab which can be
% expensive. Only stores first arriving P and S phases.
%
% INPUT
%   theFile: Name of file in which to save the TTtable. If empty, the
%            TTtable structure will not be written.
%    aModel: Name of a built-in TauP model (e.g. IASP91) or a supported TauP
%            1D model file.
%   rminmax: Vector with minimum and maximum ranges to compute travel-times
%   zminmax: Vector with minimum and maximum source depths to compute travel-times
%         N: Number of distance points to sample
%         M: Number of source depth points to sample
%
% OUTPUT
%   TTtable: Structure with the following fields
%            + phase - list of phases (currently only P and S)
%            + delta - list of depths
%            + depth - list of depths
%            + ttime - A NxMx2 array containing the travel-times at the
%                      ith-distance, jth-depth, for the kth-phase.
%            + rayP - A NxMx2 array of ray parameters
%            + inc - A NxMx2 array of ray incidince
%
% NOTES
% + The TauP toolbox has a function for creating travel-time tables that
%   is faster than making many individual calls to TauP_Time.
%

% Pre-allocate structure
TTtable.phase = {'P','S'};
TTtable.delta = linspace(rminmax(1),rminmax(2),N)';
TTtable.depth = linspace(zminmax(1),zminmax(2),M);
TTtable.ttime = zeros(N,M,2);
TTtable.rayP  = zeros(N,M,2);
TTtable.inc   = zeros(N,M,2);

tic
for ii = 1:N
    for jj = 1:M
        % Call TauP
        [ttime,aPhase,rayP,inc] = call_tauP_time(aModel,num2str(TTtable.delta(ii)),...
            '-h',num2str(TTtable.depth(jj)));
        
        % Identify first arriving P and S phase
        iP  = [];
        iS  = [];
        ind = 1;
        while isempty(iP) || isempty(iS)
            if strcmp(aPhase{ind}(end),'P') && isempty(iP)
                iP = ind;
            end
            if strcmp(aPhase{ind}(end),'S') && isempty(iS)
                iS = ind;
            end
            ind = ind + 1;
        end
        
        % Store data
        TTtable.ttime(ii,jj,1) = ttime(iP);
        TTtable.rayP(ii,jj,1)  = rayP(iP);
        TTtable.inc(ii,jj,1)   = inc(iP);
        TTtable.ttime(ii,jj,2) = ttime(iS);
        TTtable.rayP(ii,jj,2)  = rayP(iS);
        TTtable.inc(ii,jj,2)   = inc(iS);
    end
end
toc

% Save file
if ~isempty(theFile)
    save(theFile,'TTtable');
end
