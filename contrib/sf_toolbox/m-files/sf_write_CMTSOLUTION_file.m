function M = sf_write_CMTSOLUTION_file(theFile,x,y,z,t,Mw,hdr,varargin)


% In CMT catalog, r = vertical, t = south, p = east
% I copied the conversions from the strike_dip_rake_to_CMTSOLUTION.c file
% included with SPECFEM.

% Define header if not supplied
if isempty(hdr)
    hdr{1} = '2000'; % Year
    hdr{2} = '01'; % Month
    hdr{3} = '01'; % Day
    hdr{4} = '00'; % Hour
    hdr{5} = '00'; % Minute
    hdr{6} = '00'; % Second
    hdr{7} = num2str(y); % Latitude or y-coordinate
    hdr{8} = num2str(x); % Longitude or x-coordinate
    hdr{9} = num2str(z); % Depth
    hdr{10} = num2str(Mw); % Body wave magnitude
    hdr{11} = num2str(Mw); % Surface wave magnitude
    hdr{12} = 'SYNTHETIC'; % Event name
end


% Define seismic moment (from SPECFEM3D manual; in dyn*cm; 1 N = 10^5 dyn)
Mo  = 10^((3*Mw/2) + 16.1);


% Define moment tensor components
if isempty(varargin)
    % Explosive source
    Mxx = Mo/sqrt(3);
    Myy = Mo/sqrt(3);
    Mzz = Mo/sqrt(3);
    Mxy = 0;
    Mxz = 0;
    Myz = 0;
elseif length(varargin) == 3
    % Double couple source
    phi    = deg2rad(varargin{1}); % Strike CW of North
    delta  = deg2rad(varargin{2}); % Dip CW of horizontal
    lambda = deg2rad(varargin{3}); % Rake CCW from horizontal
    Mo  = 10^((3*Mw/2) + 16.1); % Seismic moment from SPECFEM3D manual (dyn*cm)
    Mxx = -Mo*(sin(delta)*cos(lambda)*sin(2*phi) + sin(2*delta)*sin(lambda)*(sin(phi)^2));
    Mxy = Mo*(sin(delta)*cos(lambda)*cos(2*phi) + 0.5*sin(2*delta)*sin(lambda)*sin(2*phi));
    Mxz = -Mo*(cos(delta)*cos(lambda)*cos(phi) + cos(2*delta)*sin(lambda)*sin(phi));
    Myy = Mo*(sin(delta)*cos(lambda)*sin(2*phi) - sin(2*delta)*sin(lambda)*(cos(phi)^2));
    Myz = -Mo*(cos(delta)*cos(lambda)*sin(phi) - cos(2*delta)*sin(lambda)*cos(phi));
    Mzz = -(Mxx + Myy); %Mo*sin(2*delta)*sin(lambda);
else
    error('Incorrect number of variable input arguments.');
end
% Interestingly, it is difficult to back out the scalar moment given the 
% moment tensor and there are many suggested ways none of which seem to
% yield above input Mo. Here are a couple of definitions for reference:
%
M = [Mxx Mxy Mxz; Mxy Myy Myz; Mxz Myz Mzz];

% (1) Here is the definition the Harvard CMT catalog uses
% Mh = sqrt(sum(M(:).^2))/sqrt(2);

% (2) Here is a simple definition from Silver and Jordan, 1982. 
% Msj = sqrt(sum(eig(M).^2)/2);


% Write the file

if isfile(theFile)
    error(['The file ',theFile,' already exists!']);
end
if ~isempty(theFile)
    fid = fopen(theFile,'w');
    
    % Header
    fprintf(fid,'%s\n',['PDE ',hdr{1},' ',hdr{2},' ',hdr{3},' ',hdr{4},' ',hdr{5},...
        ' ',hdr{6},' ',hdr{7},' ',hdr{8},' ',hdr{9},' ',hdr{10},' ',hdr{11},' ',hdr{12}]);
    % Data
    fprintf(fid,'%s\n',['event name:       ',hdr{12}]);
    fprintf(fid,'%s\n','time shift:       0.0');
    fprintf(fid,'%s %f\n','half duration:     ',t);
    fprintf(fid,'%s %f\n','latorUTM:      ',y);
    fprintf(fid,'%s %f\n','longorUTM:      ',x);
    fprintf(fid,'%s %f\n','depth:      ',abs(z));
    fprintf(fid,'%s %e\n','Mrr:      ',Mzz);
    fprintf(fid,'%s %e\n','Mtt:      ',Mxx);
    fprintf(fid,'%s %e\n','Mpp:      ',Myy);
    fprintf(fid,'%s %e\n','Mrt:      ',Mxz);
    fprintf(fid,'%s %e\n','Mrp:      ',-Myz);
    fprintf(fid,'%s %e\n','Mtp:      ',-Mxy);
    
    fclose(fid);
end

% disp(['Mrr:      ',num2str(Mzz)]);
% disp(['Mtt:      ',num2str(Mxx)]);
% disp(['Mpp:      ',num2str(Myy)]);
% disp(['Mrt:      ',num2str(Mxz)]);
% disp(['Mrp:      ',num2str(-Myz)]);
% disp(['Mtp:      ',num2str(-Mxy)]);
