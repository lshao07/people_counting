function gtrack_moduleUpdate(var, num)
global trackerProcObj;
global points;
TRACK_STATE_FREE = 0;
if trackerProcObj.gtrackHandle.activeList.count == 0
    return;
end
%tmp_list = [];
i = 1;
while i <= trackerProcObj.gtrackHandle.activeList.count
        tElem = trackerProcObj.gtrackHandle.activeList.list(i);
        uid = tElem;
        state = gtrack_unitUpdate(uid,var,num);
        if state == TRACK_STATE_FREE
            %tmp_list = [tmp_list,tElem];
            trackerProcObj.gtrackHandle.activeList = gtrack_listRemoveElement(trackerProcObj.gtrackHandle.activeList, tElem);
            trackerProcObj.gtrackHandle.targetNumCurrent = trackerProcObj.gtrackHandle.targetNumCurrent - 1;
            trackerProcObj.gtrackHandle.freeList = gtrack_listEnqueue(trackerProcObj.gtrackHandle.freeList, tElem);
        else
            i = i + 1;
        end
end

% if length(tmp_list) ~= 0
%     for j = 1:length(tmp_list)
%         trackerProcObj.gtrackHandle.activeList = gtrack_listRemoveElement(trackerProcObj.gtrackHandle.activeList, tmp_list(j));
%     end
% end
end