clear all;
close all;

datafilename = 'C:\Users\ZJ\PycharmProjects\peopleCounting\binData\pHistText_0.csv';

data = readPointCloudData(datafilename);
pointCloudSideInfoFromDSP = [];

global gObjDetObj;


% cfg.trackerDpuCfg = create_TrackerProc_Config();

gObjDetObj.dpuTrackerObj = DPU_TrackerProc_config();

outTrackerProc.numTargets = 0;
outTrackerProc.numIndices = 0;

for index = 1:length(data)

    if length(data(index).elev) == 0  
        continue;
    else       
        outTrackerProc = DPU_TrackerProc_process(gObjDetObj.dpuTrackerObj,length(data(index).elev), data(index), pointCloudSideInfoFromDSP, outTrackerProc);
        print("frame %d has target %d.", index, outTrackerProc.numTargets);
    end
end
                    