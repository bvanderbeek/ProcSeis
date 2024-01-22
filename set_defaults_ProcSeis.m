%% Default Parameters for ProcSeis

% Required paths
taup_jar   = 'contrib/TauP-2.4.5.jar'; % TauP jar-file
taup_tools = 'contrib/TauP_MatWrapper'; % TauP matlab wrappers
sf_toolbox = 'contrib/sf_toolbox'; % 

% Default Parameter Structure
% Data Panel
% MATLAB IRIS Fetch Example Data
% Defaults.DataDir      = 'EXAMPLE/IRISFETCH_MAT_DATA'; % Directory with seismic data
% Defaults.PickDir      = 'EXAMPLE/IRISFETCH_MAT_DATA/PICKS'; % Directory where picks will be saved
% Defaults.Ext          = 'mat'; % File extension on seismograms
% % SPECFEM Example Data
% Defaults.DataDir      = 'EXAMPLE/SPECFEM_DATA'; % Directory with seismic data
% Defaults.PickDir      = 'EXAMPLE/SPECFEM_DATA/PICKS'; % Directory where picks will be saved
% Defaults.Ext          = 'semv'; % File extension on seismograms
% Miniseed Example data
Defaults.DataDir      = 'EXAMPLE/MSEED_DATA'; % Directory with seismic data
Defaults.PickDir      = 'EXAMPLE/MSEED_DATA'; % Directory where picks will be saved
Defaults.Ext          = 'mseed'; % File extension on seismograms

Defaults.ChannelList  = 'BXX,BXY,BXZ'; % Channels to load in order given!
Defaults.RefModel     = 'iasp91'; % Reference 1D model for initial moveout correction
Defaults.RefPhase     = 'S'; % Reference phase for initial moveout correction
Defaults.MinTime      = -100; % Amount of data to load prior to predicted arrival
Defaults.MaxTime      = 100; % Amount of data to load following to predicted arrival
% Filter Panel
Defaults.MinFreq      = 1/40; % Minimum Butterworth filter frequency
Defaults.MaxFreq      = 1/10; % Maximum Butterworth filter frequency
Defaults.OrderFilter  = 2; % Order of Butterworth filter
% Display Panel
Defaults.ChannelName  = {'E','N','Z'}'; % Channel names
Defaults.ChannelIndex = 1; % Channel index to plot
Defaults.Gain         = 1e6; % Scale factor for seismograms
Defaults.TraceIndices = 1:30; % Trace indices to load
% Phase Panel
Defaults.WindowLength = 0; % Cross-correlation window length
Defaults.TaperLength  = 0; % Hanning taper length to add to ends of window
Defaults.PreBuffer    = 0; % Extend window this many seconds earlier
Defaults.PostBuffer   = 0; % Extend window this many seconds later
