function [GDT,DTF,DTS] = GenSwaveDelay(handles)
% 1) Measure delays in PAZ coordinates
% 2) Run this

% Must be in PAZ coordinate system
xyz = handles.MenuCoord.String{handles.MenuCoord.Value};
tf_fsz = false;
if ~strcmp(xyz,'PAZ')
    if strcmp(xyz,'FSZ')
        tf_fsz = true;
        disp('Using current coordinates...');
    else
        error('Must run from PAZ coordinates!');
    end
end

% Define parameters
t       = handles.SeisDat.t - handles.picks.t_stack;
ichn    = handles.MenuChan.Value;
wlength = str2double(handles.DTxtWLen.String);
tbuff   = [str2double(handles.DTxtPreBuff.String),str2double(handles.DTxtPostBuff.String)];
ttaper  = str2double(handles.DTxtWTap.String);
maxlag = round(handles.SeisDat.fs*wlength/2);

% Define amplitude scaling
if handles.ChkNormAmp.Value
    % Get RMS trace amplitudes in phase window
    A = get_rms_amp(t,squeeze(handles.SeisDat.seis(:,ichn,:)),handles.picks.DT,wlength);
else
    A = ones(handles.SeisDat.ntrace,1);
end
A = A(~handles.picks.tf_bad);
A = repmat(A,1,handles.SeisDat.nsamp);


% Create stacked signal components
u1 = squeeze(handles.SeisDat.seis(~handles.picks.tf_bad,1,:));
u1 = stack_seis(t,u1./A,handles.picks.DT(~handles.picks.tf_bad));
u2 = squeeze(handles.SeisDat.seis(~handles.picks.tf_bad,2,:));
u2 = stack_seis(t,u2./A,handles.picks.DT(~handles.picks.tf_bad));
if handles.SeisDat.nchan == 3
    u3   = squeeze(handles.SeisDat.seis(~handles.picks.tf_bad,3,:));
    u3   = stack_seis(t,u3./A,handles.picks.DT(~handles.picks.tf_bad));
    Spaz = cat(1,u1,u2,u3);
else
    Spaz = cat(1,u1,u2);
end
A    = mean(mean(A,1),2);
Spaz = A*Spaz;

% Compute cross-correlations with respect to reference stacked waveform
GDT = zeros(handles.SeisDat.ntrace,1);
DTF = zeros(handles.SeisDat.ntrace,1);
DTS = zeros(handles.SeisDat.ntrace,1);
for itr = 1:handles.SeisDat.ntrace
    
    if tf_fsz
        S   = Spaz;
        FSZ = squeeze(handles.SeisDat.seis(itr,:,:));
    else
        % (1) Rotate back to geographic (ENZ) coordinates
        % Prior rotation matrix
        R = reshape(handles.SeisDat.Rmat(itr,:),handles.SeisDat.nchan,handles.SeisDat.nchan);
        % Reverse prior rotation for ith-trace
        FSZ = (R')*squeeze(handles.SeisDat.seis(itr,:,:));
        % Define reference waveform for ith-trace
        S = (R')*Spaz;
        
        % (2) Rotate to Fast-Slow coordinates
        % Index station
        [~,ista] = ismember(handles.SeisDat.station(itr),handles.station.id);
        % Get stored fast axis orientation
        alpha = handles.station.rotation.FSZ(ista);
        % Rotate fast direction to channel 1 from ENZ
        R = [cos(alpha),-sin(alpha), 0;...
            sin(alpha), cos(alpha), 0;...
            0,          0,          1];
        R = R(1:handles.SeisDat.nchan,1:handles.SeisDat.nchan);
        % Apply rotation (undoing the rotation hence the transpose)
        FSZ = (R')*FSZ;
        % Apply rotation to reference waveform
        S = (R')*S;
    end
    
    % (3) Window waveforms
    w   = get_window(handles.SeisDat.nsamp,handles.SeisDat.fs,ttaper,t(1),...
        -tbuff(1),wlength + sum(tbuff));
    wfs = get_window(handles.SeisDat.nsamp,handles.SeisDat.fs,ttaper,t(1),...
        handles.picks.DT(itr)-tbuff(1),wlength + sum(tbuff));
    S   = S.*repmat(w(:)',size(S,1),1);
    FSZ = FSZ.*repmat(wfs(:)',size(FSZ,1),1);
    
    % (4) Compute cross-correlations
    [rf,lags] = xcorr(FSZ(1,:),S(1,:),maxlag,'coeff');
    % Note that the cross-correlation for the slow waveform is reveresed 
    % with respect to fast waveform
    rs        = xcorr(S(2,:),FSZ(2,:),maxlag,'coeff');
    lags      = lags./handles.SeisDat.fs;
    
    % (5) Define generalized delay time
    r        = rf + rs;
    keep     = (r == max(r));
    GDT(itr) = mean(lags(keep));
    
    keep     = (rf == max(rf));
    DTF(itr) = mean(lags(keep));
    keep     = (rs == max(rs));
    DTS(itr) = mean(-lags(keep));
    
%     keep     = (rf == max(rf));
%     dtf      = mean(lags(keep));
%     keep     = (rs == max(rs));
%     dts      = mean(lags(keep));
%     if (dtf > -dts) && (abs(GDT(itr)) > 5)
%         rf       = xcorr(-FSZ(1,:),S(1,:),maxlag,'coeff');
%         r        = rf + rs;
%         keep     = (r == max(r));
%         GDT(itr) = mean(lags(keep));
%     end
end
