function value = sf_read_input(theFile,keyword,delimeter,tf_num)
% SF_READ_INPUT: Simple function to read values from SPECFEM3D style input
% file. Basically just a system call to grep. This function searches the
% input file name for keyword using grep and returns all text found on the
% line containing the keyword occuring after the specified delimeter but 
% excluding any text occuring after a # comment character.
%

% Fix spaces in file name as this causes a problem for grep
theFile = strsplit(theFile,' ');
if length(theFile) > 1
    for ii = 1:(length(theFile)-1)
        theFile{ii} = [theFile{ii},'\ '];
    end
    theFile = cat(2,theFile{:});
else
    theFile = theFile{1};
end

[stat,value] = system(['grep ''^',keyword,''' ',theFile]);
if stat == 0
    % Process value
    in = min(strfind(value,delimeter)) + length(delimeter);
    im = min(strfind(value,'#')) - 1;
    if isempty(im)
        im = length(value);
    end
    % Remove white space
    value = deblank(value(in:im)); % trailing white space
    value = deblank(value(end:-1:1)); % leading white space
    value = value(end:-1:1);
    % Convert to number?
    if tf_num
        value = str2double(value);
        if isnan(value)
            error('Problem interpreting value from file.');
        end
    end
else
    error('Call to grep returned with an error.');
end