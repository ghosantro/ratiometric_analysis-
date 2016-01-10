%Return the data in a .dat in a three-dimensional array, where
%the third coordinate is the frame number.
%
%opread(file) returns the whole .dat file as a uint16 matrix.
%opread(file,precision) returns the matrix as array of class 'precision'.
%opread(file, start) returns everything after the frame start as a matrix.
%opread(file, start, finish) returns everything between the frame start and
%   the frame finish as a matrix.
%opread(file,start,finish,[x y]) returns time series from frame start to
%   finish for pixel (x,y).

% Modifications:
% 2006-03-02 initialize data matrix with X333 = zeros(header(4), header(3), ...);
%            (change of x and y) SL
%            introduced ncol/nrow taking into account binning SL
% 2006-06-07 added single pixel read capability
%            enabled by passing a vector [x,y] as the fourth argument to
%            the function.
%            was not verified with binned images.  -AS
% 2007-03-13 array v is initialized as a uint16 rather than default double
% 2007-04-28 added precision option, with uint16 as default
% 2007-05-15 can specify pixel without specifying start and finish

function [v1,v2,ncol,nrow] = Dual_imaging_opread(file,dual_image, start, finish)
% Dual_image is a 0 or 1 index 
% 0 = single imaging 
% 1 = dual imaging 

if (nargin<1)
    error('opread requires at least a filename.');
end

header = opheader(file);
precision = 'uint16';
if start ==0
    start=1;
end
if finish ==0
    finish=header{2};
end 

xbin = header{5};
ybin = header{6};

% ncol = header{3}/xbin;
% nrow = header{4}/ybin;

ncol = header{3}
nrow = header{4}

% nparam = 0;
% 
% for i = 1:nargin-1
%     if ischar(varargin{i}) precision = varargin{i}; % read precision
%     elseif length(varargin{i})>1
%         pixel = varargin{i}; % pixel [x,y]
%     elseif nparam == 0
%         start=varargin{i};
%         nparam = nparam+1;
%     elseif nparam == 1
%         finish=varargin{i};
%         nparam = nparam+1;
%     end
% end

%flip_image_choice = varargin{end};

machineformat = 'ieee-le'; % 'ieee-be'

clear fid
fid = fopen(file, 'rb', machineformat);
version = char(fread(fid, 1, 'char'));%Determine the version we're using.


if dual_image==0
    v1=zeros(ncol, nrow, finish-start+1,precision); %initialize the data matrix.
    v2 = []; 
    fseek(fid, 0, -1);

    if version=='a' || version=='c' || version=='d'
        fseek(fid, 1024, -1);
        fseek(fid, 2*nrow*ncol*(start-1), 0);
        for loop = 1:(finish-start+1)         
             X3 = fread(fid, ncol*nrow, 'uint16', 'l'); %Each block of data has size width*height         
             X33 = reshape(X3, ncol, nrow); %Turn it into a square matrix.
    %          if flip_image_choice==1
    %             elseif flip_image_choice==2
    %                  X33=fliplr(X33);
    %             elseif flip_image_choice==3
    %                  X33=flipud(X33);
    %             elseif flip_image_choice==4
    %                  X33=rot90(X33,2);
    %             elseif flip_image_choice==5
    %                  X33=rot90(X33,3);
    %             elseif flip_image_choice==6
    %                  X33=rot90(X33,1);
    %             end
              v1(:,:,loop)=X33; %Write the result as the next frame in the data matrix.
         end
        fclose(fid);
    %     else % else offset to pixel's position ...
    %         offset = 2*( (pixel(1)-1)*ncol + (ncol-pixel(2))); %data is arrayed in reverse order, so must subtract from max indices
    %         fseek(fid,offset,0);
    %         for loop = 1:(finish-start+1) % ... and read one pixel from each frame
    %             v(loop) = fread(fid, 1, 'uint16', 2*(ncol*nrow-1), 'l');
    %         end
    %         fclose(fid);
    %     end
        
    else
        disp ('Not a recognized version');
    end 
elseif dual_image==1
    v1=zeros(ncol, nrow, finish-start+1,precision); %initialize the data matrix.
    v2 =zeros(ncol, nrow, finish-start+1,precision); 
    fseek(fid, 0, -1);

    if version=='a' || version=='c' || version=='d'
        fseek(fid, 1024, -1);
        fseek(fid, 2*nrow*ncol*(start-1), 0);
        for loop = 1:(finish-start+1)         
             X3 = fread(fid, ncol*nrow, 'uint16', 'l'); %Each block of data has size width*height         
             X33 = reshape(X3, ncol, nrow);
             X33_show = double(X33);%Turn it into a square matrix.
            
    %          if flip_image_choice==1
    %             elseif flip_image_choice==2
    %                  X33=fliplr(X33);
    %             elseif flip_image_choice==3
    %                  X33=flipud(X33);
    %             elseif flip_image_choice==4
    %                  X33=rot90(X33,2);
    %             elseif flip_image_choice==5
    %                  X33=rot90(X33,3);
    %             elseif flip_image_choice==6
    %                  X33=rot90(X33,1);
    %             end
              if loop/2~=floor(loop/2)
                  v1(:,:,loop)=X33; %Write the result as the next frame in the data matrix.
              else 
                  v2(:,:,loop) = X33;
                  
              end 
         end
        fclose(fid);
        
        %post-processing the two images interpolate the points in between the real  
        if floor((finish-start+1)/2)~=(finish-start+1)/2
            v1(:,:,2:2:end)=(v1(:,:,1:2:end-2)+v1(:,:,3:2:end))./2;
            v2(:,:,1) = v2(:,:,2); 
            v2(:,:,end) = v2(:,:,end-1); 
            v2(:,:,3:2:end-1)= (v2(:,:,2:2:end-2)+v2(:,:,4:2:end))./2; 
        elseif floor((finish-start+1)/2)==(finish-start+1)/2
            v2(:,:,1) = v2(:,:,2);
            v2(:,:,3:2:end) = (v2(:,:,2:2:end-2)+v2(:,:,4:2:end))./2; 
            v1(:,:,end) = v1(:,:,end-1); 
            v1(:,:,2:2:end-1)= (v1(:,:,1:2:end-3)+v1(:,:,3:2:end))./2; 
                  
        end
        
    %     else % else offset to pixel's position ...
    %         offset = 2*( (pixel(1)-1)*ncol + (ncol-pixel(2))); %data is arrayed in reverse order, so must subtract from max indices
    %         fseek(fid,offset,0);
    %         for loop = 1:(finish-start+1) % ... and read one pixel from each frame
    %             v(loop) = fread(fid, 1, 'uint16', 2*(ncol*nrow-1), 'l');
    %         end
    %         fclose(fid);
    %     end
        
    else
        disp ('Not a recognized version');
    end 
end 

return;