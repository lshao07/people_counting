function radarData = readPointCloudData(filename)
radarData = [];
fileID = fopen(filename);
headerline = strtrim(fgetl(fileID));

while ~feof(fileID)
    frameInfo = split(strtrim(fgetl(fileID)), {',', ', '});
    pointCloud = split(strtrim(fgetl(fileID)), {',', ', '});
    trackList = split(strtrim(fgetl(fileID)), {',', ', '});
    
    frameInfo = frameInfo(1:end-1);
    pointCloud = pointCloud(1:end-1);
    trackList = trackList(1:end-1);
    
    frameStruct = struct();
    
    %add frame number to frameStruct
    frameStruct.frameNum = frameInfo(17);
    
    %add point cloud to frameStruct
    frameStruct.elev = [];
    frameStruct.azim = [];
    frameStruct.doppler = [];
    frameStruct.range = [];
    frameStruct.snr = [];
    for n = 1:length(pointCloud)/5 - 1
        frameStruct.elev(end+1) = str2double(pointCloud(1+n*5));
        frameStruct.azim(end+1) = str2double(pointCloud(2+n*5));
        frameStruct.doppler(end+1) = str2double(pointCloud(3+n*5));
        frameStruct.range(end+1) = str2double(pointCloud(4+n*5));
        frameStruct.snr(end+1) = str2double(pointCloud(5+n*5));      
    end
    
    %add track list to frameStruct
    frameStruct.TID = [];
    frameStruct.x = [];
    frameStruct.y = [];
    frameStruct.z = [];
    frameStruct.vx = [];
    frameStruct.vy = [];
    frameStruct.vz = [];
    frameStruct.ax = [];
    frameStruct.ay = [];
    frameStruct.az = [];
    for n = 1:length(trackList)/10 - 1
        frameStruct.TID(end+1) = str2double(trackList(1+n*10));
        frameStruct.x(end+1) = str2double(trackList(2+n*10));
        frameStruct.y(end+1) = str2double(trackList(3+n*10));
        frameStruct.z(end+1) = str2double(trackList(4+n*10));
        frameStruct.vx(end+1) = str2double(trackList(5+n*10));
        frameStruct.vy(end+1) = str2double(trackList(6+n*10));
        frameStruct.vz(end+1) = str2double(trackList(7+n*10));
        frameStruct.ax(end+1) = str2double(trackList(8+n*10));
        frameStruct.ay(end+1) = str2double(trackList(9+n*10));
        frameStruct.az(end+1) = str2double(trackList(10+n*10));    
    end
    
    %add single frameStruct to totalData
    radarData = cat(1, radarData, frameStruct);
   
end

fclose(fileID);

end