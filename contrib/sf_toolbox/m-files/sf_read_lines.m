function params = sf_read_lines(theFile,varargin)
% SF_READ_LINES: Reads an input file line by line and returns each line in
% a cell array of strings.

if isempty(varargin)
    nadd = 0;
else
    nadd = varargin{1};
end

[stat,nmax] = system(['wc -l ',theFile]);
if stat == 0
    nmax = strsplit(strtrim(nmax),' ');
    nmax = str2double(nmax{1});
    params = cell(nmax,1);
end
if ~(stat == 0) || isnan(nmax)
    error(['Problem reading number of lines in ',theFile]);
end

% Read file line by line
fid = fopen(theFile);
n   = 0; % Line number counter
m   = 0; % Non-commented line counter
while n < nmax
    % Read past lines that are commented
    char1 = '#';
    while strcmp(char1,'#')
        nline = fgetl(fid);
        ichar = find(~isspace(nline),1,'first');
        char1 = nline(ichar);
        n = n + 1;
    end
    m = m + 1; % Counts number of non-commented lines
    params{m} = nline;
    % params{m} = str2num(nline); %#ok <nline may be an array>
end
fclose(fid);
params = params(1:m);

