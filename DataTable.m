%% Data Table Fields

Station           = []; % Station Name
Event             = []; % Event ID
Delta             = []; % Station-Event Distance (deg.)
BAZ               = []; % Station-Event Backazimuth (deg.)
Phase             = []; % Seismic Phase
Channel           = []; % Channel Analysed
Orientation       = []; % Azimuth and elevation of Channel
Filter            = []; % Filter type
Order             = []; % Order of filter
ZeroPhase         = []; % Logical, true is zero-phase filter used
CornerFrequencies = []; % Corner frequencies of filter (Hz)
ReferenceModel    = []; % Reference velocity model
DataType          = []; % Data type measured (e.g., Arrival Time, Splitting Intensity, Splitting Parameters)
PredictedValue    = []; % Predicted value from reference model
MeasuredValue     = []; % Measured Value
Uncertainty       = []; % Uncertainty in measurement
Method            = []; % Measurement method (e.g., MCC)
Window            = []; % Window length used in measurement (s)
