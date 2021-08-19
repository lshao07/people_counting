function gtrack_moduleAllocate(num)
global trackerProcObj;
global points;
% trackerProcObj.gtrackHandle = trackerProcObj.gtrackHandle;
GTRACK_ID_POINT_NOT_ASSOCIATED = 255;
isSnrThresholdPassed = false;
for n = 1:num
    if trackerProcObj.gtrackHandle.bestIndex(n) == GTRACK_ID_POINT_NOT_ASSOCIATED
        trackerProcObj.gtrackHandle.allocIndex  = n;
        allocnum = 1;
        allocSNR = points(n).snr;
        
        tElemFree = trackerProcObj.gtrackHandle.freeList.list(1);
        mCenter.vector = points(n);
        mSum_vector = points(n);
        for k = n+1:num %循环范围感觉有点问题 
            if k ~= n
                if trackerProcObj.gtrackHandle.bestIndex(k) == GTRACK_ID_POINT_NOT_ASSOCIATED
                    mCurrent.vector = points(k);
                    mCurrent.vector.doppler = gtrack_unrollRadialVelocity(trackerProcObj.gtrackHandle.params.maxRadialVelocity, mCenter.vector.doppler, mCurrent.vector.doppler);
                    if abs(mCurrent.vector.doppler - mCenter.vector.doppler) < trackerProcObj.gtrackHandle.params.allocationParams.maxVelThre
                        distance = gtrack_calcDistance(mCenter.vector,mCurrent.vector);
                        if sqrt(distance) < trackerProcObj.gtrackHandle.params.allocationParams.maxDistanceThre
                            trackerProcObj.gtrackHandle.allocIndex = cat(2,trackerProcObj.gtrackHandle.allocIndex,k);
                            allocnum = allocnum + 1;
                            allocSNR = allocSNR +  points(k).snr;
                            %add mCenter to mSum
                            mSum_vector.azimuth = mSum_vector.azimuth + mCurrent.vector.azimuth;
                            mSum_vector.doppler = mSum_vector.doppler + mCurrent.vector.doppler;
                            mSum_vector.elevation = mSum_vector.elevation + mCurrent.vector.elevation;
                            mSum_vector.range = mSum_vector.range + mCurrent.vector.range;
                            mSum_vector.snr = mSum_vector.snr + mCurrent.vector.snr;
                            
                            mCenter.vector.azimuth = mSum_vector.azimuth/allocnum;
                            mCenter.vector.doppler = mSum_vector.doppler/allocnum;
                            mCenter.vector.elevation = mSum_vector.elevation/allocnum;
                            mCenter.vector.range = mSum_vector.range/allocnum;
                            mCenter.vector.snr = mSum_vector.snr/allocnum;
                        end
                    end
                end
            end
        end
        
        if allocnum > trackerProcObj.gtrackHandle.params.allocationParams.pointsThre && abs(mCenter.vector.doppler) > trackerProcObj.gtrackHandle.params.allocationParams.velocityThre
            isBehind = false;
            if trackerProcObj.gtrackHandle.activeList.count ~= 0
                for i = 1:trackerProcObj.gtrackHandle.activeList.count
                    tElem = trackerProcObj.gtrackHandle.activeList.list(i);
                    uid = tElem;
                    hs = trackerProcObj.gtrackHandle.hTrack(uid).H_s.array;
                    if(gtrack_isPointBehindTarget(mCenter.vector, hs))  
                        isBehind = true;
                        break;
                    end
                end
            end
            if (isBehind)
                isSnrThresholdPassed = allocSNR > trackerProcObj.gtrackHandle.params.allocationParams.snrThreObscured;
            else
                isSnrThresholdPassed = allocSNR > trackerProcObj.gtrackHandle.params.allocationParams.snrThre;
            end
            
            if isSnrThresholdPassed
                for m = 1:allocnum
                    trackerProcObj.gtrackHandle.bestIndex(trackerProcObj.gtrackHandle.allocIndex(m)) = tElemFree;
                end
                
                %allocate new tracker
                trackerProcObj.gtrackHandle.targetNumTotal = trackerProcObj.gtrackHandle.targetNumTotal + 1;
                trackerProcObj.gtrackHandle.targetNumCurrent = trackerProcObj.gtrackHandle.targetNumCurrent + 1;
                %Dequeue
                if gtrack_isListEmpty(trackerProcObj.gtrackHandle.freeList)
                    return
                end
                [tElemFree,trackerProcObj.gtrackHandle.freeList] = gtrack_listDequeue(trackerProcObj.gtrackHandle.freeList);
                uid = tElemFree;
                gtrack_unitStart(uid,trackerProcObj.gtrackHandle.heartBeat, trackerProcObj.gtrackHandle.targetNumTotal, mCenter.vector);
                
                %Enqueue
                trackerProcObj.gtrackHandle.activeList = gtrack_listEnqueue(trackerProcObj.gtrackHandle.activeList, tElemFree);
            end
        end
    end
end
% trackerProcObj.gtrackHandle = trackerProcObj.gtrackHandle;

end