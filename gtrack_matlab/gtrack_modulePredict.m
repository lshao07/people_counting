function gtrack_modulePredict()
    global trackerProcObj;
    if trackerProcObj.gtrackHandle.activeList.count == 0
        return;
    end
    for i = 1:trackerProcObj.gtrackHandle.activeList.count
        tElem = trackerProcObj.gtrackHandle.activeList.list(i);
        uid = tElem;
        gtrack_unitPredict(uid);
    end
end