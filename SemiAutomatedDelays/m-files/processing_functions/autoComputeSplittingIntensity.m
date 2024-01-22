function F = autoComputeSplittingIntensity(evtid,ATT,TT1D,event,station,dataDir,refModel)

% Define analysis parameters
sampFreq     = 10; %40; % Sampling frequency (Hz)
timeWindow   = 120; %600; % Length of seismograms (s)
corners      = [1/40, 1/10]; %[1/33, 1/12]; % Corner frequencies for Butterworth filter (Hz)
order        = 2; %3; % Order of Butterworth filter
tf_zerophase = true; % Use zero-phase filter?
% Splitting Intensity Analysis Parameters
tstart  = 0; % Window start time (s)
tlength = 15; % Window length (s)
inorm   = 3; %1;
% Window limits for polarisation recomputation
tf_pol = false;
pmin   = 0;
pmax   = 15;
% Arrival identification parameters
CHN = 'PZ'; % Original channel picked
FLT = [0, 1/40, 1/10, 2, 0]; %[0,1/33,1/12,3,0]; % Original filter used
WND = 8; %10; % Only original picks should have a window parameter of 0 s
tf_specfem = true;

% Get phase
aPhase = unique(ATT.phase(ATT.event == evtid));
if length(aPhase) == 1
    aPhase = aPhase{1};
elseif isempty(aPhase)
    error('Event has no pre-existing picks!');
else
    if evtid == 697
        warning('Assuming S for event 697.');
        aPhase = 'S';
    else
        error('Multiple phases recorded by this event!');
    end
end

if tf_specfem
    sfEvent = sf_read_events(dataDir);
    sfStation = sf_read_stations(dataDir);
    data = sf_load_bin_seis(sfStation,sfEvent,evtid,{'BXX','BXY','BXZ'},dataDir,'semv',timeWindow*[-0.5, 0.5],...
        'filtopts',[corners(1),corners(2),order,tf_zerophase],'tf_rotate',true,...
        'tf_moveout',true,'ttap',[],'aPhase',aPhase,'aModel',refModel);
    % Modify data structure
    data.t = data.t(:)';
    data = rmfield(data,'event');
    data.event.id         = evtid;
    data.event.originTime = datenum(event.originDate(event.id == evtid)); % DATENUM
    data.event.latitude   = event.latitude(event.id == evtid);
    data.event.longitude  = event.longitude(event.id == evtid);
    data.event.depth      = event.depth(event.id == evtid);
    if data.ntrace == length(station.name)
        data.nsta = data.ntrace;
        data = rmfield(data,'ntrace');
    else
        error('Missing traces!');
    end
    [~,lis] = ismember(data.station,station.name);
    data.network = station.network(lis);
    data.refTime = data.event.originTime + (data.tt1D./(60*60*24));
    data.seis = permute(data.seis,[1,3,2]);
else
    data = loadSeis_commonEvent(dataDir,evtid,event,station,sampFreq,...
        timeWindow,refModel,aPhase,corners,order,tf_zerophase);
end

%% Subset Seismograms: Only use those with picks
keep = (ATT.event == evtid) & strcmp(ATT.channel,CHN)...
    & ismember(ATT.filter,FLT,'rows') & (ATT.window == WND);

igood = ismember(data.station,unique(ATT.station(keep)));
data.network = data.network(igood);
data.station = data.station(igood);
data.nsta    = sum(igood);
data.delta   = data.delta(igood);
data.baz     = data.baz(igood);
data.refTime = data.refTime(igood);
data.seis    = data.seis(igood,:,:);

%% Compute Delay Times

% Compute mean station delays
% Relative delay times
keep = strcmp(ATT.channel,CHN) & ismember(ATT.filter,FLT,'rows') & (ATT.window == WND);
sub_eid = unique(ATT.event(keep));
sub_sid = unique(ATT.station(keep));
[~,ievt] = ismember(ATT.event(keep),sub_eid);
[~,ista] = ismember(ATT.station(keep),sub_sid);
% Absolute delay times
ADT = 24*60*60*(datenum(ATT.arrivalTime(keep)) - datenum(ATT.originTime(keep))) - TT1D(keep);
% Mean event delays
MED = accumarray(ievt,ADT,[],@mean,NaN);
% Relative delat times
RDT = ADT - MED(ievt);
% Relative Mean station delays
MSD = accumarray(ista,RDT,[],@mean,NaN);
% Polarisation
pazm = ATT.pol(keep,1);
pelv = ATT.pol(keep,2);

% Associate Delays to Seismograms
[~,kevt] = ismember(data.event.id,sub_eid);
[iin,ksta] = ismember(data.station,sub_sid);
[~,karr] = ismember([kevt*ones(data.nsta,1),ksta],[ievt,ista],'rows');

% Check for unique polarisation
pazm = unique(pazm(karr));
pelv = unique(pelv(karr));
if (length(pazm) > 1) || (length(pelv) > 1)
    error('Non-unique polarisation!');
end

if ~(sum(iin) == length(iin))
    error('Unexpected behaviour');
end

%% Recompute Polarisation

if tf_pol
    % Polarisation window
    pwin = (data.t >= pmin) & (data.t <= pmax);
    % Align TRZ Traces
    TRZ = data.seis;
    TRZ = align_traces(data.t,TRZ,MED(kevt) + RDT(karr));
    % Scale by 3-component amplitude in window
    TRZ = scale_traces(pwin,TRZ,[]);
    % Identify polarity reversals
    T     = mean(squeeze(TRZ(:,:,1)),1);
    iflip = TRZ(:,pwin,1)*(T(:,pwin)')./sqrt(sum(TRZ(:,pwin,1).^2,2).*sum(T(:,pwin).^2,2));
    iflip = iflip < -0.1; % Allow a little bit of anti-correlation
    % Correct suspected polarity reversals; assumes all channels are effected!
    TRZ(iflip,:,:) = -TRZ(iflip,:,:);
    % Stack
    TRZ = squeeze(mean(TRZ,1));
    % Estimate Polarisation Azimuth and elevation
    [pazm,pelv,D,V] = compute_polarization(squeeze(TRZ(pwin,1)),squeeze(TRZ(pwin,2)),squeeze(TRZ(pwin,3)));
    % Find rotation angle about x-axis to minimize P-wave energy
    R = roty(pelv)*rotz(-pazm);
    V = R*V;
    a = 90 - atan2d(V(3,1),V(2,1));
    % Final rotation matrix
    R = rotx(a)*R;
    
%     % Using directly the permuted Eigenvectors was giving polarity
%     % reversals that modified the splitting intensity measurements.
%     % Sort as (1) principle S-wave, (2) orthogonal S-wave, (3) P-wave
%     % Principle S-wave is that with the largest magnitude eigenvalue
%     i1 = 3;
%     % Identify P-wave as component with most energy on z-component
%     i3 = find(abs(V(3,1:2)) == max(abs(V(3,1:2))),1,'first');
%     % Orthogonal S-wave is remaining component
%     if i3 == 1
%         i2 = 2;
%     elseif i3 == 2
%         i2 = 1;
%     end
%     % Rotation matrix that maps from TRZ to PNL coordinates
%     R = (V(:,[i1,i2,i3]))';
else
    R = roty(pelv)*rotz(-pazm);
end

% % Plot Check
% PNL = (R*(TRZ'))'; % Sesmograms rotated to principle coordinates
% H = figure;
% subplot(1,2,1); hold on;
% plot(data.t,TRZ(:,1),'-b','linewidth',2);
% plot(data.t,TRZ(:,2),'-r','linewidth',2);
% plot(data.t,TRZ(:,3),'-g','linewidth',2);
% plot([pmin,pmin],[-1,1],'--k','linewidth',2);
% plot([pmax,pmax],[-1,1],'--k','linewidth',2);
% box on; grid on; axis tight;
% xlim([pmin-tlength,pmax+tlength]);
% subplot(1,2,2); hold on;
% plot(data.t,PNL(:,1),'-b','linewidth',2);
% plot(data.t,PNL(:,2),'-r','linewidth',2);
% plot(data.t,PNL(:,3),'-g','linewidth',2);
% plot([pmin,pmin],[-1,1],'--k','linewidth',2);
% plot([pmax,pmax],[-1,1],'--k','linewidth',2);
% box on; grid on; axis tight;
% xlim([pmin-tlength,pmax+tlength]);
% H.Position(3) = 2*H.Position(3);

%% Rotate Seismograms

% Define principle S-wave
Spnl = zeros(data.nsta,data.nsamp,3);
for ii = 1:data.nsta
    % Rotate seismograms
    si = squeeze(data.seis(ii,:,:))';
    si = R*si;
    % Store principle component
    Spnl(ii,:,:) = si';
end

%% Splitting Intensity
% A bit confused on how to apply trace amplitude scaling. Particularly
% important when constructing stacked trace.
ttaper = 0; % HARD-CODED NO TAPER ON WINDOW

% Get window for aligned traces
w    = get_window(data.nsamp,data.fs,ttaper,data.t(1),tstart,tlength)';
iwin = w > 0;

% Aligned and scale traces. Scaling is inconsequential for splitting
% intensity calculations as it is applied to all components.
S0     = align_traces(data.t,Spnl,MED(kevt) + RDT(karr));
[S0,A] = scale_traces(iwin,S0,[]); % SCALE BY 3-COMPONENT AMPLITUDE
Spnl   = Spnl./repmat(A,1,data.nsamp,3);

% Identify polarity reversals
s1 = mean(squeeze(S0(:,:,1)),1);
% Identify polarity reversals
iflip = S0(:,iwin,1)*(s1(:,iwin)')./sqrt(sum(S0(:,iwin,1).^2,2).*sum(s1(:,iwin).^2,2));
iflip = iflip < -0.1; % Allow a little bit of anti-correlation
% Correct suspected polarity reversals
S0(iflip,:,:) = -S0(iflip,:,:);
Spnl(iflip,:,:) = -Spnl(iflip,:,:);

% Derivative of stacked principle component trace
s1   = mean(squeeze(S0(:,:,1)),1);
s2   = mean(squeeze(S0(:,:,2)),1);
dsdt = 0.5*data.fs*(s1(3:end) - s1(1:end-2));
dsdt = cat(2,data.fs*(s1(2) - s1(1)),dsdt);
dsdt = cat(2,dsdt,-data.fs*(s1(end-1) - s1(end)));

% Derivative of individual traces
U1   = squeeze(S0(:,:,1));
U2   = squeeze(S0(:,:,2));
dudt = 0.5*data.fs*(U1(:,3:end) - U1(:,1:end-2));
dudt = cat(2,data.fs*(U1(:,2) - U1(:,1)),dudt);
dudt = cat(2,dudt,-data.fs*(U1(:,end-1) - U1(:,end)));

% % Plot check alignment
% H = figure;
% subplot(1,2,1); hold on;
% plot(data.t,squeeze(S0(:,:,1)));
% plot(data.t,s1,'-k','linewidth',2);
% plot([tstart*ones(2,1), (tstart+tlength)*ones(2,1)],...
%     max(abs(s1(:)))*[-1,-1; 1,1],'-g','linewidth',2);
% box on; grid on; axis tight; xlim([-max(1./corners), tlength + max(1./corners)]);
% subplot(1,2,2); hold on;
% plot(data.t,squeeze(S0(:,:,2)));
% plot(data.t,dsdt,'-k','linewidth',2);
% plot([tstart*ones(2,1), (tstart+tlength)*ones(2,1)],...
%     max(abs(s1(:)))*[-1,-1; 1,1],'-g','linewidth',2);
% box on; grid on; axis tight; xlim([-max(1./corners), tlength + max(1./corners)]);
% H.Position(3) = 2*H.Position(3);

% Compute trace magnitudes. I think this should be done before
% application of window because if the window contains a taper it could
% bias the magnitude calculation.
% Principle component magnitudes in window
as1 = sqrt(sum(s1(iwin).^2));
au1 = sqrt(sum(U1(:,iwin).^2,2));
% Total magnitudes in window
As = sqrt(sum((s1(iwin).^2) + (s2(iwin).^2)));
Au = sqrt(sum((U1(:,iwin).^2) + (U2(:,iwin).^2),2));

% Apply window to principle component derivatives
dsdt = iwin.*dsdt;
dudt = repmat(iwin,data.nsta,1).*dudt;

% Shift traces by mean station delay + mean event delay
tq = repmat(data.t,data.nsta,1) - MED(kevt) - repmat(MSD(ksta),1,data.nsamp);
dsdt = interp1(data.t,dsdt,tq,'linear',0);
for jj = 1:data.nsta
    dudt(jj,:) = interp1(data.t,dudt(jj,:),tq(jj,:),'linear',0);
end
% Need to interpolate window for computing errors
iwin = interp1(data.t,w,tq,'linear',0);
iwin = iwin > 0;

% Option 2: How to normalize trace amplitudes?
% Scaling reference trace to account for amplitude loss from splitting.
% This is necessary because the reference trace should be that which is
% uneffected by splitting.
switch inorm
    case 1
        % No normalization. Nothing to do.
    case 2
        % Normalize to stacked trace three-component magnitude
        dsdt = (As/as1)*dsdt;
        dudt = repmat(As./au1,1,data.nsamp).*dudt;
    case 3
        % Normalize to individual trace three component magnitude
        dsdt = repmat(Au./as1,1,data.nsamp).*dsdt;
        dudt = repmat(Au./au1,1,data.nsamp).*dudt;
end

% Compute splitting intensity
SI  = sum(squeeze(Spnl(:,:,2)).*dsdt,2)./sum(dsdt.*dsdt,2);
SIj = sum(squeeze(Spnl(:,:,2)).*dudt,2)./sum(dudt.*dudt,2);

% Estimated errors
% Only consider piece of waveform in window
T = squeeze(Spnl(:,:,2));
T = T.*iwin;
SI  = cat(2,SI,sum((T - dsdt.*repmat(SI,1,data.nsamp)).^2,2)./sum(dsdt.*dsdt,2));
SIj = cat(2,SIj,sum((T - dudt.*repmat(SIj,1,data.nsamp)).^2,2)./sum(dudt.*dudt,2));

% Results Structure
F.eventid = data.event.id*ones(data.nsta,1);
F.station = data.station;
F.dlt     = data.delta;
F.baz     = data.baz;
F.paz     = pazm*ones(data.nsta,1);
F.SI      = SI;
F.SIj     = SIj;

