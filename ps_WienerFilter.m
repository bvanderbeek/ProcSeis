function handles = ps_WienerFilter(handles)

% Get Parameters
ichn          = handles.MenuChan.Value;
s             = eval(handles.DTxtGain.String);
tf_norm       = handles.ChkNormAmp.Value;
t_stack       = handles.picks.t_stack;
dt            = handles.picks.DT;
fs            = handles.SeisDat.fs;
nsamp         = handles.SeisDat.nsamp;
winlen        = str2double(handles.DTxtWLen.String);
ttaper        = str2double(handles.DTxtWTap.String);
tbuff         = [str2double(handles.DTxtPreBuff.String),str2double(handles.DTxtPostBuff.String)];
xmax          = str2double(handles.DTxtXmax.String);

% Phase arrival time vector
t = handles.SeisDat.t - t_stack;

% Define amplitude scaling
A = ones(1,handles.SeisDat.ntrace);
if tf_norm
    % Get RMS trace amplitudes in phase window
    if winlen <= 0
        dtnw = xmax;
    else
        dtnw = winlen;
    end
    
    % Window each trace and compute RMS
    for ii = 1:handles.SeisDat.ntrace
        % lwin  = (t >= 0) & (t <= (dtnw - dt(ii))); % Was wrong?
        lwin  = ((t - dt(ii)) >= 0) & ((t - dt(ii)) <= dtnw);
        A(ii) = rms(squeeze(handles.SeisDat.seis(ii,ichn,lwin)));
    end
end
% Interpolate traces to align phases
W = squeeze(handles.SeisDat.seis(:,ichn,:))';
U = interp2(1:handles.SeisDat.ntrace,t(:),W,repmat(1:handles.SeisDat.ntrace,handles.SeisDat.nsamp,1),...
    repmat(t,1,handles.SeisDat.ntrace) + repmat(dt(:)',handles.SeisDat.nsamp,1),'linear',0);
% Remove bad traces
U(:,handles.picks.tf_bad) = 0;
% Apply amplitude scaling
U = U./repmat(A,handles.SeisDat.nsamp,1);
% Stack
U = mean(A(~handles.picks.tf_bad))*sum(U,2)./sum(~handles.picks.tf_bad);

% Apply tapered window to stacked waveform
win = get_window(nsamp,fs,ttaper,t(1),-tbuff(1),winlen + sum(tbuff));
U   = win.*U;

% Also window waveforms?
win = interp1(t(:)',win(:)',repmat(t,1,handles.SeisDat.ntrace) - repmat(dt(:)',handles.SeisDat.nsamp,1),'linear',0);
W   = win.*W;

% Frequency-domain Wiener Filter
U  = fft(U,[],1);
W  = fft(W,[],1);
dw = repmat(0.001*max(abs(W.*W),[],1),nsamp,1);
F  = abs(W.*U./max(W.*W,dw));
W  = real(ifft(F.*W,[],1));
% U = real(ifft(U,[],1));

% Update seismograms
handles.SeisDat.seis(:,ichn,:) = permute(repmat(W',1,1,1),[1,3,2]);
