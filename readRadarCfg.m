function [radarCfg] = readRadarCfg(cfgFileName)
    radarCfg = struct();
    fileID = fopen(cfgFileName);
    while ~feof(fileID)
        tline = strtrim(fgetl(fileID));
        cfg = split(string(tline), {'   ', '  ', ' '});
        if length(cfg) > 1
           name = cfg(1);
           value = str2double(cfg(2:end));
           radarCfg = setfield(radarCfg,name,value);
        else
%            radarCfg = cat(1, radarCfg, cfg); 
        end       
    end
    fclose(fileID);
end