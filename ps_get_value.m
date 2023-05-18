function varargout = ps_get_value(handles,varargin)

nout      = length(varargin);
varargout = cell(1,nout);

for iout = 1:nout
    switch varargin{iout}
        % DATA PANEL
        case 'DataDir'
            varargout{iout} = handles.DTxtDataDir.String;
        case 'PickDir'
            varargout{iout} = handles.DTxtPickDir.String;
        case 'ChannelList'
            varargout{iout} = strsplit(handles.DTxtChan.String,',');
        case 'Ext'
            varargout{iout} = handles.DTxtExt.String;
        case 'RefModel'
            varargout{iout} = handles.DTxtRefMod.String;
        case 'RefPhase'
            varargout{iout} = handles.DTxtPhase.String;
        case 'MinTime'
            varargout{iout} = str2double(handles.DTxtMinTime.String);
        case 'MaxTime'
            varargout{iout} = str2double(handles.DTxtMaxTime.String);
        % EVENT PANEL
        case 'EventIndex'
            varargout{iout} = str2double(handles.DTxtEind.String);
        % FILTER PANEL
        case 'MinFreq'
            varargout{iout} = str2double(handles.DTxtMinFreq.String);
        case 'MaxFreq'
            varargout{iout} = str2double(handles.DTxtMaxFreq.String);
        case 'OrderFilter'
            varargout{iout} = str2double(handles.DTxtOrder.String);
        % DISPLAY PANEL
        case 'XMin'
            varargout{iout} = str2double(handles.DTxtXmin.String);
        case 'XMax'
            varargout{iout} = str2double(handles.DTxtXmax.String);
        case 'Gain'
            varargout{iout} = str2double(handles.DTxtGain.String);
        case 'TraceIndices'
            varargout{iout} = eval(handles.DTxtNTrace.String);
        case 'ChannelName'
            varargout{iout} = handles.MenuChan.String{handles.MenuChan.Value};
        case 'ChannelIndex'
            varargout{iout} = handles.MenuChan.Value;
        case 'CoordName'
            varargout{iout} = handles.MenuCoord.String{handles.MenuCoord.Value};
        case 'CoordIndex'
            varargout{iout} = handles.MenuCoord.Value;
        case 'SortName'
            varargout{iout} = handles.MenuSort.String{handles.MenuSort.Value};
        case 'SortIndex'
            varargout{iout} = handles.MenuSort.Value;
        % PHASE PANEL
        case 'WindowLength'
            varargout{iout} = str2double(handles.DTxtWLen.String);
        case 'TaperLength'
            varargout{iout} = str2double(handles.DTxtWTap.String);
        case 'PreBuffer'
            varargout{iout} = str2double(handles.DTxtPreBuff.String);
        case 'PostBuffer'
            varargout{iout} = str2double(handles.DTxtPostBuff.String);
        otherwise
            error(['Unknown variable ',varargin{iout},'!']);
    end
end