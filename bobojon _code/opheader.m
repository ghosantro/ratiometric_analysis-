% OPHEADER(FILE) returns header information from surface map file
%   HEADER = OPHEADER(FILE) returns
%       header{1} - the date and time the file was created, as a string,
%       header{2} - the number of frames,
%       header{3} - the width,
%       header{4} - the height,
%    all as doubles.
%
%   Version 'd' files also return
%       header{5} - xbin,
%       header{6} - ybin.
%
%   Written by Martin , Cornell University, 2005.

% Modifications:
% 2006-03-02 SL Introduced version 'd' to include xbin - header{5}
%                                                 ybin - header{6}

function v = opheader(file)

machineformat = 'ieee-le';
permission = 'rb';

fid = fopen(file,permission,machineformat);
version = char(fread(fid,1, 'char'));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HEADER VERSION A
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (version=='a')
    datetime1 = fread(fid, 20, 'char');
    datetime = (cast(transpose(datetime1), 'char'));
    nframes= fread(fid, 1, 'double');
    height= fread(fid, 1, 'double');
    width= fread(fid, 1, 'double');
    flag=0;
    comment='';
    fseek(fid, 1, 0);
    while ~flag
        x=fread(fid, 1, 'char');
        if x==3 flag=1;
        else comment = [comment x];
        end
    end
    v = {datetime, nframes, height, width, comment};
    return;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HEADER VERSION C
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (version=='c')
    datetime1 = fread(fid, 17, 'char');
    datetime = (cast(transpose(datetime1), 'char'));
    nframes= fread(fid, 1, 'int');
    height= fread(fid, 1, 'int');
    width= fread(fid, 1, 'int');
    flag=0;
    comment='';
    fseek(fid, 1, 0);
    while ~flag
        x=fread(fid, 1, 'char');
        if x==3 flag=1;
        else comment = [comment x];
        end
    end
    v = {datetime, nframes, height, width, comment};
    return;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HEADER VERSION D
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (version=='d')
    
    datetime1 = fread(fid, 17, 'char');  %early files are 17, later are 24
    % if date string begins with an alphabet character ('A' = 65) then the
    % date string is 24 characters long and thus 7 more need to be read.
    if datetime1(1) >= 65
        datetime1 = [datetime1; fread(fid, 7, 'char')];
    end
    datetime = (cast(transpose(datetime1), 'char'));
    nframes = fread(fid, 1, 'int');
    height = fread(fid, 1, 'int');
    width = fread(fid, 1, 'int');
    xbin = fread(fid, 1, 'int');
    ybin = fread(fid, 1, 'int');
    mark1 = fread(fid, 1, 'char');
    mark3 = fread(fid, 1, 'char');
    flag=0;
    comment='';
    fseek(fid, 1, 0);
    flag = 1;
    while ~flag
        x=char(fread(fid, 1, 'char'));
        if x==3 flag=1;
        else comment = [comment x];
        end
    end
    v = {datetime, nframes, height, width, xbin, ybin, comment};
    return;
    
else
    disp 'Not a recognized version';
end

fclose(fid);

return;