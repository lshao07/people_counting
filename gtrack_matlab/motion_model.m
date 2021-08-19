clear all;
close all;

addpath 'D:\matlab_workspace\people_count';
datafilename = 'D:\matlab_workspace\people_count\pHistText_0.csv';

global trackerProcObj;
global outTrackerProc;

data = readPointCloudData(datafilename);
pointCloudSideInfoFromDSP = [];


% cfg.trackerDpuCfg = create_TrackerProc_Config();

trackerProcObj = DPU_TrackerProc_config();


outTrackerProc.numTargets = 0;
outTrackerProc.numIndices = 0;

for index = 1:length(data)

    if length(data(index).elev) == 0  
        continue;
    else       
        targetList = DPU_TrackerProc_process(length(data(index).elev), data(index), pointCloudSideInfoFromDSP);
        if outTrackerProc.numTargets == length(data(index).TID)
            
        else
            fprintf("frame %d has %d targets, algorithm computes %d targets.\n", index, length(data(index).TID), outTrackerProc.numTargets);            
            if length(data(index).TID) ~= 0
                fprintf("data_target_tid is: ")
                disp(data(index).TID);
            end
            if isfield(targetList, 'tid')
                fprintf("algo_target_tid is: ")
                disp([targetList.tid]);
            end
        end
        
    end
end
                    