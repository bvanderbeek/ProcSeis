function handles = ps_plot_stack(handles)
% Toggle axes flag
ps_switch_axis_stat(handles.AxStat1,2);

% Get Plotting Parameters
ichn          = handles.MenuChan.Value;
s             = eval(handles.DTxtGain.String);
tf_norm       = handles.ChkNormAmp.Value;
t_stack       = handles.picks.t_stack;
dt            = handles.picks.DT;
winlen        = str2double(handles.DTxtWLen.String);
xmin          = str2double(handles.DTxtXmin.String);
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
u = squeeze(handles.SeisDat.seis(:,ichn,:))';
u = interp2(1:handles.SeisDat.ntrace,t(:),u,repmat(1:handles.SeisDat.ntrace,handles.SeisDat.nsamp,1),...
    repmat(t,1,handles.SeisDat.ntrace) + repmat(dt(:)',handles.SeisDat.nsamp,1),'linear',0);
% Remove bad traces
u(:,handles.picks.tf_bad) = 0;
% Apply amplitude scaling
u = u./repmat(A,handles.SeisDat.nsamp,1);
% Stack
u = s*mean(A(~handles.picks.tf_bad))*sum(u,2)./sum(~handles.picks.tf_bad);

cla(handles.AxStack);
plot(handles.AxStack,t,u,'-k','linewidth',2);

% Plot phase window
if winlen > 0
    hold(handles.AxStack,'on');
    
    % Additional window parameters
    dt_buff = [str2double(handles.DTxtPreBuff.String),str2double(handles.DTxtPostBuff.String)];
    dt_tap  = str2double(handles.DTxtWTap.String);
    
    % Main window
    xbox = [0,winlen,winlen,0,0];
    ybox = [min(u),min(u),max(u),max(u),min(u)];
    fill(handles.AxStack,xbox,ybox,'g','EdgeColor','None','FaceAlpha',0.25);
    % Pre-arrival buffer
    if dt_buff(1) > 0
        xbox = [-dt_buff(1),0,0,-dt_buff(1),-dt_buff(1)];
        fill(handles.AxStack,xbox,ybox,'y','EdgeColor','None','FaceAlpha',0.25);
    end
    % Post-arrival buffer
    if dt_buff(2) > 0
        xbox = winlen + [0,dt_buff(2),dt_buff(2),0,0];
        fill(handles.AxStack,xbox,ybox,'y','EdgeColor','None','FaceAlpha',0.25);
    end
    % Taper
    if dt_tap > 0
        xbox = -dt_buff(1) + [-dt_tap,0,0,-dt_tap,-dt_tap];
        plot(handles.AxStack,xbox,ybox,'--k','linewidth',2);
        xbox = winlen + dt_buff(2) + [0,dt_tap,dt_tap,0,0];
        plot(handles.AxStack,xbox,ybox,'--k','linewidth',2);
    end
end

% Update Axes
axis(handles.AxStack,'tight');
handles.AxStack.XLim  = [xmin,xmax];
handles.AxStack.Box   = 'on';
handles.AxStack.XGrid = 'on';
handles.AxStack.YGrid = 'on';

% Toggle axes flag
ps_switch_axis_stat(handles.AxStat1,0);