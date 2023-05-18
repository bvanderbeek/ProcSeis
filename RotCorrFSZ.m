function handles = RotCorrFSZ(handles)

% Define parameters
t       = handles.SeisDat.t - handles.picks.t_stack;
wlength = str2double(handles.DTxtWLen.String);
tbuff   = [str2double(handles.DTxtPreBuff.String),str2double(handles.DTxtPostBuff.String)];
ttaper  = str2double(handles.DTxtWTap.String);

% Identify split based on rotation correlation
handles.picks.dt_split  = zeros(handles.SeisDat.ntrace,1);
handles.picks.fast_azim = zeros(handles.SeisDat.ntrace,1);
for itr = 1:handles.SeisDat.ntrace
    % (1) Define seismogram for ith trace
    S = squeeze(handles.SeisDat.seis(itr,:,:));
    
    % (2) Window waveforms
    w = get_window(handles.SeisDat.nsamp,handles.SeisDat.fs,ttaper,t(1),...
        handles.picks.DT(itr)-tbuff(1),wlength + sum(tbuff));
    S = S.*repmat(w(:)',size(S,1),1);
    
    % (3) Find fast-axis direction
    [alpha,dts] = splitRC(S(1:2,:),handles.SeisDat.fs,0.4*wlength,181,false);
    alpha = alpha(1) - sign(dts(1))*(pi/2);
    dts   = dts(1);
    % (4) Rotate to Fast-Slow coordinates
    % Rotate fast direction to channel 1 from ENZ
    R = [cos(alpha),-sin(alpha), 0;...
        sin(alpha), cos(alpha), 0;...
        0,          0,          1];
    R = R(1:handles.SeisDat.nchan,1:handles.SeisDat.nchan);
    % Apply rotation (undoing the rotation hence the transpose)
    % S = (R')*S;
    handles.SeisDat.seis(itr,:,:) = (R')*squeeze(handles.SeisDat.seis(itr,:,:));
%     if dts > 0
%         handles.SeisDat.seis(itr,1,:) = -handles.SeisDat.seis(itr,1,:);
%     end
    
    % (5) Store Results
    % Splitting parameters
    handles.picks.dt_split(itr)  = dts;
    handles.picks.fast_azim(itr) = alpha;
    % Rotation matrix
    R0 = reshape(handles.SeisDat.Rmat(itr,:),handles.SeisDat.nchan,handles.SeisDat.nchan);
    R  = (R')*R0;
    handles.SeisDat.Rmat(itr,:) = R(:)';
    
end
