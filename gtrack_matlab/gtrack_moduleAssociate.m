function gtrack_moduleAssociate(num)
    global trackerProcObj;
    global points;
    if trackerProcObj.gtrackHandle.activeList.count == 0
        return;
    end
    for i = 1:trackerProcObj.gtrackHandle.activeList.count
        tElem = trackerProcObj.gtrackHandle.activeList.list(i);
        uid = tElem;
        gtrack_unitScore(uid,num);
    end
end