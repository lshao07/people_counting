function [tNum,targetDescr] = gtrack_moduleReport(tNum,targetDescr)
    global trackerProcObj;
    if trackerProcObj.gtrackHandle.activeList.count == 0
        return;
    end
    for i = 1:trackerProcObj.gtrackHandle.activeList.count
        tElem = trackerProcObj.gtrackHandle.activeList.list(i);
        uid = tElem;
        target = gtrack_unitReport(uid);
        tNum = tNum + 1;
        targetDescr = cat(1, targetDescr, target);
    end
end