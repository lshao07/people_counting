function gtrack_unitEvent(uid, num)
global trackerProcObj;

TRACK_STATE_DETECTION = 2;
TRACK_STATE_ACTIVE = 3;
TRACK_STATE_FREE = 0;

if trackerProcObj.gtrackHandle.hTrack(uid).state == TRACK_STATE_DETECTION
    if num > trackerProcObj.gtrackHandle.hTrack(uid).allocationParams.pointsThre
        % Hit Event 
		trackerProcObj.gtrackHandle.hTrack(uid).detect2freeCount = 0;
		trackerProcObj.gtrackHandle.hTrack(uid).detect2activeCount = trackerProcObj.gtrackHandle.hTrack(uid).detect2activeCount + 1;
        if(trackerProcObj.gtrackHandle.hTrack(uid).detect2activeCount > trackerProcObj.gtrackHandle.hTrack(uid).stateParams.det2actThre)
            trackerProcObj.gtrackHandle.hTrack(uid).state = TRACK_STATE_ACTIVE;
        end
    else
        if num == 0
            % Miss 
			trackerProcObj.gtrackHandle.hTrack(uid).detect2freeCount = trackerProcObj.gtrackHandle.hTrack(uid).detect2freeCount + 1;
            if(trackerProcObj.gtrackHandle.hTrack(uid).detect2activeCount > 0)
                trackerProcObj.gtrackHandle.hTrack(uid).detect2activeCount = trackerProcObj.gtrackHandle.hTrack(uid).detect2activeCount - 1;
            end
            if(trackerProcObj.gtrackHandle.hTrack(uid).detect2freeCount > trackerProcObj.gtrackHandle.hTrack(uid).stateParams.det2freeThre)
                trackerProcObj.gtrackHandle.hTrack(uid).state = TRACK_STATE_FREE;
            end
        else
            trackerProcObj.gtrackHandle.hTrack(uid).detect2freeCount = 0;
        end
    end
elseif trackerProcObj.gtrackHandle.hTrack(uid).state == TRACK_STATE_ACTIVE
    if num
        % Hit Event 
		trackerProcObj.gtrackHandle.hTrack(uid).active2freeCount = 0;
    else
        % Miss
        trackerProcObj.gtrackHandle.hTrack(uid).active2freeCount = trackerProcObj.gtrackHandle.hTrack(uid).active2freeCount + 1;
        
        % Set threshold based on whether target is static, or is in exit zone */
        % If target is static and inside the static box */
        if(trackerProcObj.gtrackHandle.hTrack(uid).sceneryParams.numStaticBoxes)
            thre = trackerProcObj.gtrackHandle.hTrack(uid).stateParams.exit2freeThre;
            for numBoxes = 1:trackerProcObj.gtrackHandle.hTrack(uid).sceneryParams.numStaticBoxes
                if(gtrack_isPointInsideBox(trackerProcObj.gtrackHandle.hTrack(uid).S_hat, trackerProcObj.gtrackHandle.hTrack(uid).sceneryParams.staticBox))
                    if(trackerProcObj.gtrackHandle.hTrack(uid).isTargetStatic)
                        thre = trackerProcObj.gtrackHandle.hTrack(uid).stateParams.static2freeThre;
                    else
                        thre = trackerProcObj.gtrackHandle.hTrack(uid).stateParams.active2freeThre;
                    end
                    break;
                end
            end
        else
            % Normal moving target */
            thre = trackerProcObj.gtrackHandle.hTrack(uid).stateParams.active2freeThre;
        end
        
        % Threshold can not be more than lifetime of the target */
        if(thre > trackerProcObj.gtrackHandle.hTrack(uid).heartBeatCount)
            thre = trackerProcObj.gtrackHandle.hTrack(uid).heartBeatCount;
        end
        if(trackerProcObj.gtrackHandle.hTrack(uid).active2freeCount > thre) 
            trackerProcObj.gtrackHandle.hTrack(uid).state = TRACK_STATE_FREE;
        end
    end
else
end

end