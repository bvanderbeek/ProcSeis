function handles = ps_set_value(handles,Name,Value)

% Allow for character array and cell array of character vectors
if ~iscell(Name)
    Name  = {Name};
    Value = {Value};
end

nout = length(Name);
for iout = 1:nout
    switch Name{iout}
        % DATA PANEL
        case 'DataDir'
            handles.DTxtDataDir.String = Value{iout};
        case 'PickDir'
            handles.DTxtPickDir.String = Value{iout};
        case 'ChannelList'
            handles.DTxtChan.String = Value{iout};
        case 'Ext'
            handles.DTxtExt.String = Value{iout};
        case 'RefModel'
            handles.DTxtRefMod.String = Value{iout};
        case 'RefPhase'
            handles.DTxtPhase.String = Value{iout};
        case 'MinTime'
            handles.DTxtMinTime.String = num2str(Value{iout});
        case 'MaxTime'
            handles.DTxtMaxTime.String = num2str(Value{iout});
        % EVENT PANEL
        case 'EventIndex'
            handles.DTxtEind.String = num2str(Value{iout});
        % FILTER PANEL
        case 'MinFreq'
            handles.DTxtMinFreq.String = num2str(Value{iout});
        case 'MaxFreq'
            handles.DTxtMaxFreq.String = num2str(Value{iout});
        case 'OrderFilter'
            handles.DTxtOrder.String = num2str(Value{iout});
        % DISPLAY PANEL
        case 'XMin'
            handles.DTxtXmin.String = num2str(Value{iout});
        case 'XMax'
            handles.DTxtXmax.String = num2str(Value{iout});
        case 'Gain'
            handles.DTxtGain.String = num2str(Value{iout});
        case 'TraceIndices'
            if mean(diff(Value{iout})) == 1
                handles.DTxtNTrace.String = [num2str(Value{iout}(1)),':',num2str(Value{iout}(end))];
            else
                handles.DTxtNTrace.String = num2str(Value{iout});
            end
        case 'ChannelName'
            handles.MenuChan.String = Value{iout};
        case 'ChannelIndex'
            handles.MenuChan.Value = Value{iout};
        case 'CoordName'
            handles.MenuCoord.String = Value{iout};
        case 'CoordIndex'
            handles.MenuCoord.Value = Value{iout};
        case 'SortName'
            handles.MenuSort.String = Value{iout};
        case 'SortIndex'
            handles.MenuSort.Value = Value{iout};
        % PHASE PANEL
        case 'WindowLength'
            handles.DTxtWLen.String = num2str(Value{iout});
        case 'TaperLength'
            handles.DTxtWTap.String = num2str(Value{iout});
        case 'PreBuffer'
            handles.DTxtPreBuff.String = num2str(Value{iout});
        case 'PostBuffer'
            handles.DTxtPostBuff.String = num2str(Value{iout});
        otherwise
            error(['Unknown variable ',Name{iout},'!']);
    end
end
