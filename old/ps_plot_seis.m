function Ax = ps_plot_seis(Ax,data,trace_indices,sort_var,ichan,s,dt,tf_fill,tf_norm)

limt    = [];

% Sort index
if isempty(sort_var)
    sort_var = 'traceid';
    data.traceid = (1:data.ntrace)';
end
[~,ind] = sort(data.(sort_var));

% Make figures
trace_indices = fliplr(trace_indices(:)');
hold(Ax,'on');
for jj = trace_indices % Trace indices correspond to sorted order
    % Arrival-adjusted time vector
    t = data.t - dt(ind(jj));
    
    % Define scaled and normalized (optional) trace
    if ~isempty(limt) && tf_norm
        f = squeeze(data.seis(ind(jj),ichan,:));
        f = 2*f./rms(f((t >= limt(1)) & (t <= limt(2))));
    else
        f = s*squeeze(data.seis(ind(jj),ichan,:));
    end
    
    % Plot traces with option to add waveform filling
    if tf_fill
        % Filled waveforms
        n = length(f);
        f = cat(1,f,zeros(n,1));
        t = cat(1,t,flipud(t));
        fill(Ax,t,jj + max(f,0),'b','EdgeColor','None','FaceAlpha',0.5);
        fill(Ax,t,jj + min(f,0),'r','EdgeColor','None','FaceAlpha',0.5);
        plot(Ax,t(1:n),jj + f(1:n),'-k');
    else
        % Line waveforms
        plot(Ax,t,jj + f,'-k');
    end
end
% Define title
ipage = ceil(max(trace_indices)/length(trace_indices));
npage = 1 + floor(data.ntrace/length(trace_indices));
Ax.Title.String = ['Channel ',data.channel{ichan},'----Event ',num2str(data.event),...
    '----Traces ',num2str(min(trace_indices)),' to ',num2str(max(trace_indices)),...
    ' (page ',num2str(ipage),' of ',num2str(npage),')'];

% Set up labels
Ax.YTick         = sort(trace_indices);
YTickLabel       = data.(sort_var)(ind(trace_indices));
Ax.YTickLabel    = cellstr(num2str(YTickLabel(:)));
Ax.YLabel.String = sort_var;
axis(Ax,'tight');
