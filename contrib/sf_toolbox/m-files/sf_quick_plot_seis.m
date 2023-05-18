function H = sf_quick_plot_seis(ichan,data,s,sort_var,dt,n2plot,varargin)

limt    = [];
tf_norm = false;
tf_fill = false;
if length(varargin) == 1
    limt = varargin{1};
elseif length(varargin) == 2
    limt    = varargin{1};
    tf_norm = varargin{2};
elseif length(varargin) == 3
    limt    = varargin{1};
    tf_norm = varargin{2};
    tf_fill = varargin{3};
end

% Number of figures needed
N       = unique([1:n2plot:data.ntrace,data.ntrace]);
% Sort index
if isempty(sort_var)
    sort_var = 'traceid';
    data.traceid = (1:data.ntrace)';
end
[~,ind] = sort(data.(sort_var));

% Make figures
H = cell(length(N)-1,1);
for ii = 1:(length(N)-1)
    H{ii} = figure('Units','normalized','Position',[0,0,0.4,1]);
    for jj = N(ii):N(ii+1)
        hold on;
%         if ~isempty(limt) && tf_norm
%             t = data.t - dt(ind(jj));
%             f = squeeze(data.seis(ind(jj),ichan,:));
%             limt = varargin{1};
%             f = 2*f./rms(f((t >= limt(1)) & (t <= limt(2))));
%             
%             plot(t,jj + f,'-k');
%         else
%             if tf_fill
%                 % Filled waveforms
%                 f = [s*squeeze(data.seis(ind(jj),ichan,:));0;0];
%                 t = data.t - dt(ind(jj));
%                 t = cat(1,t,t(end),t(1));
%                 fill(t,jj + max(f,0),'b','EdgeColor','None','FaceAlpha',0.5);
%                 fill(t,jj + min(f,0),'r','EdgeColor','None','FaceAlpha',0.5);
%                 plot(t,jj + f,'-k');
%             else
%                 plot(data.t - dt(ind(jj)),jj + s*squeeze(data.seis(ind(jj),ichan,:)),'-k');
%             end
%         end
        t = data.t - dt(ind(jj));
        if ~isempty(limt) && tf_norm
            f = squeeze(data.seis(ind(jj),ichan,:));
            limt = varargin{1};
            f = 2*f./rms(f((t >= limt(1)) & (t <= limt(2))));
        else
            f = s*squeeze(data.seis(ind(jj),ichan,:));
        end
        if tf_fill
            % Filled waveforms
            f = cat(1,f,0,0);
            t = cat(1,t,t(end),t(1));
            fill(t,jj + max(f,0),'b','EdgeColor','None','FaceAlpha',0.5);
            fill(t,jj + min(f,0),'r','EdgeColor','None','FaceAlpha',0.5);
            plot(t(1:end-2),jj + f(1:end-2),'-k');
        else
            plot(t,jj + f,'-k');
        end
    end
    box on;
    grid on;
    xlabel('time');
    title(['Channel ',num2str(ichan),'; Traces ',num2str(N(ii)),' to ',num2str(N(ii+1))]);
    
    % Set up labels
    mlab = 5;
    H{ii}.Children.YTick = N(ii):mlab:N(ii+1);
    YLabels = data.(sort_var)(ind(N(ii):N(ii+1)));
    YLabels = YLabels(1:mlab:end);
    H{ii}.Children.YTickLabel = cellstr(num2str(YLabels(:)));
    ylabel(sort_var);
    
    if ~isempty(limt)
        xlim(limt);
    end
end

