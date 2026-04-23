function OUT = BarVideo(Maps,startpos,barspeed,angle,videolength)
%Maps Maps of the rhabdomere
%startpos [x y width height] the start pos of the bar
%barspeed speed of the bar left to right
%videolegth number of frames in video
xdim = size(Maps,2);
ydim = size(Maps,1);
nmaps =size(Maps,3);
OUT =zeros(videolength,nmaps);
minx = min(startpos(:,1));
maxx = max(startpos(:,1)+startpos(:,3));
miny = min(startpos(:,2));
maxy = max(startpos(:,2)+startpos(:,4));
frame = zeros(maxy-miny+1,maxx-minx+1);
for i =1:size(startpos,1)
    % wrec = ones(startpos(i,4),startpos(i,3)); %square
    %ellipse
    wrec = zeros(startpos(i,4),startpos(i,3));
    hhalf =startpos(i,4)/2;
    whalf =startpos(i,3)/2;
    [Xc,Yc] = meshgrid((0:(startpos(i,3)-1))+0.5-whalf,(0:(startpos(i,4)-1))+0.5-hhalf);
    halfgrid = 0.5/(startpos(i,4)+startpos(i,3))*2;
    wrec((Xc/whalf).^2+(Yc/hhalf).^2<=(1-halfgrid)) =1;
   
    wrect = imtranslate(wrec,[startpos(i,1)-minx startpos(i,2)-miny], 'OutputView','full');
    frame(1:size(wrect,1),1:size(wrect,2)) =frame(1:size(wrect,1),1:size(wrect,2)) + wrect;
end
for i =1:videolength
     image = zeros(ydim,xdim);
     xpos =round(minx + (i-1)*barspeed);
    
    framet =imtranslate(frame,[xpos miny], 'OutputView','full');
    if(xpos > 0)
        xpos =0;
    end
    framec =imcrop(framet, [-xpos 0 xdim ydim]);
    image(1:size(framec,1),1:size(framec,2)) =framec;
    
    image =imrotate(image,angle,'nearest','crop');
    for j = 1:nmaps
        OUT(i,j) =sum(sum(image.*squeeze(Maps(:,:,j))));
    end
end
end

