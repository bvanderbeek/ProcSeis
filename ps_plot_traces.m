function handles = ps_plot_traces(handles)
ps_switch_axis_stat(handles.AxStat2,2);

% Get Plotting Parameters
ichan         = handles.MenuChan.Value;
s             = eval(handles.DTxtGain.String);
sort_var      = handles.MenuSort.String{handles.MenuSort.Value};
trace_indices = eval(handles.DTxtNTrace.String);
tf_fill       = handles.ChkFillSeis.Value;
tf_norm       = handles.ChkNormAmp.Value;
t_stack       = handles.picks.t_stack;
dt            = handles.picks.DT;
winlen        = str2double(handles.DTxtWLen.String);
xmin          = str2double(handles.DTxtXmin.String);
xmax          = str2double(handles.DTxtXmax.String);

% Window parameters
if winlen > 0
    dt_buff = [str2double(handles.DTxtPreBuff.String),str2double(handles.DTxtPostBuff.String)];
    % dt_tap  = str2double(handles.DTxtWTap.String);
    xbox    = [-dt_buff(1),winlen+dt_buff(2),winlen+dt_buff(2),-dt_buff(1),-dt_buff(1)];
    ybox    = 0.4*[-1,-1,1,1,-1];
end

% Define Normalization Window
if tf_norm
    if winlen <= 0
        dtnw = xmax;
    else
        dtnw = winlen;
    end
%     % Normalization factor based on rms-amplitude in stack
%     lwin = (handles.AxStack.Children(end).XData > 0) & (handles.AxStack.Children(end).XData < dtnw);
%     arms = rms(handles.AxStack.Children(end).YData(lwin))/s;
    
    % Normalization factor based on rms-amplitude in all windows
    lwin = ((handles.SeisDat.t(:)') >= (t_stack + dt(:))) &...
        ((handles.SeisDat.t(:)') <= (t_stack + dt(:) + dtnw));
    arms = squeeze(handles.SeisDat.seis(:,ichan,:));
    arms = rms(arms(lwin));
end

% Define Trace Sorting
if strcmp(sort_var,'traceid')
    sort_var = (1:handles.SeisDat.ntrace)';
else
    sort_var = handles.SeisDat.(sort_var);
end
% Arrange traces in ascending order of 'sort_var'
[~,ind] = sort(sort_var,'ascend');

% Need to reverse order of trace_indices to correctly label traces because
% the first sorted trace is plotted at the top of the figure window.
trace_indices = fliplr(trace_indices(:)');

% Plot Traces
cla(handles.AxTraces,'reset');
hold(handles.AxTraces,'on');
for jj = trace_indices % Trace indices correspond to sorted order
    % Arrival-adjusted time vector
    t = handles.SeisDat.t - t_stack - dt(ind(jj));
    
    % Define scaled and normalized (optional) trace
    if tf_norm
        lwin = (t >= 0) & (t <= dtnw);
        f = squeeze(handles.SeisDat.seis(ind(jj),ichan,:));
        f = s*arms*f./rms(f(lwin));
    else
        f = s*squeeze(handles.SeisDat.seis(ind(jj),ichan,:));
    end
    
    % The DC shift for plotting multiple traces defined such that the first
    % trace in the sort is plotted at the top of the figure window.
    df = trace_indices(1) - jj + trace_indices(end);
    
    % Plot traces with option to add waveform filling
    if tf_fill && ~handles.picks.tf_bad(ind(jj))
        % Define waveform closed shape
        n = length(f);
        f = cat(1,f,zeros(n,1));
        t = cat(1,t,flipud(t));
        % Plot filled waveforms
        fill(handles.AxTraces,t,max(f,0) + df,'b','EdgeColor','None','FaceAlpha',0.5);
        fill(handles.AxTraces,t,min(f,0) + df,'r','EdgeColor','None','FaceAlpha',0.5);
        plot(handles.AxTraces,t(1:n),f(1:n) + df,'-k');
    elseif ~tf_fill && ~handles.picks.tf_bad(ind(jj))
        % Line waveforms
        plot(handles.AxTraces,t,f + df,'-k');
    else
        % Bad/Disregarded trace
        plot(handles.AxTraces,t,f + df,'-r','linewidth',2);
    end
    
    % Plot phase window
    if (winlen > 0) && ~handles.picks.tf_bad(ind(jj))
        fill(handles.AxTraces,xbox,ybox + df,[0.5,0.5,0.5],'EdgeColor','None','FaceAlpha',0.25);
    end
    
end
% Define title
ipage = ceil(max(trace_indices)/length(trace_indices));
npage = 1 + floor(handles.SeisDat.ntrace/length(trace_indices));
handles.AxTraces.Title.String = ['Channel ',handles.MenuChan.String{ichan},'----Event ',num2str(handles.SeisDat.event),...
    '----Traces ',num2str(min(trace_indices)),' to ',num2str(max(trace_indices)),...
    ' (page ',num2str(ipage),' of ',num2str(npage),')'];

% Define Trace Labels
YTickLabel = sort_var(ind(trace_indices));

% Update Axes
% Trace labels
handles.AxTraces.YTick         = sort(trace_indices);
handles.AxTraces.YTickLabel    = cellstr(num2str(YTickLabel(:)));
handles.AxTraces.YLabel.String = handles.MenuSort.String{handles.MenuSort.Value};
axis(handles.AxTraces,'tight');
% Station Names
yyaxis(handles.AxTraces,'right');
handles.AxTraces.YTick      = (handles.AxTraces.YAxis(1).TickValues - handles.AxTraces.YAxis(1).Limits(1))./diff(handles.AxTraces.YAxis(1).Limits);
handles.AxTraces.YTickLabel = handles.SeisDat.station(ind(trace_indices));
handles.AxTraces.YColor     = [0,0.4470,0.7410];
yyaxis(handles.AxTraces,'left');
% Plot box
handles.AxTraces.XLim  = [xmin,xmax];
handles.AxTraces.Box   = 'on';
handles.AxTraces.XGrid = 'on';
handles.AxTraces.YGrid = 'on';

ps_switch_axis_stat(handles.AxStat2,0);