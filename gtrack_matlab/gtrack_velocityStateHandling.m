function uVec = gtrack_velocityStateHandling(uid,center)
global trackerProcObj;

inst = trackerProcObj.gtrackHandle.hTrack(uid);
rvIn = center(4);
if inst.velocityHandling == 0
    uVec = inst.rangeRate;
    inst.velocityHandling = 1;
elseif inst.velocityHandling == 1
    instanteneousRangeRate = (center(1) - inst.allocationRange)/(inst.heartBeatCount-inst.allocationTime)*inst.dt;
    
    inst.rangeRate = inst.unrollingParams.alpha * inst.rangeRate + (1-inst.unrollingParams.alpha) * instanteneousRangeRate;
    uVec = gtrack_unrollRadialVelocity(inst.maxRadialVelocity, inst.rangeRate, rvIn);
    
    rrError = (instanteneousRangeRate - inst.rangeRate)/inst.rangeRate;
    
    if abs(rrError) < inst.unrollingParams.confidence
        inst.velocityHandling = 2;
    end
elseif inst.velocityHandling == 2
    instanteneousRangeRate = (center(1) - inst.allocationRange)/(inst.heartBeatCount-inst.allocationTime)*inst.dt;
    
    inst.rangeRate = inst.unrollingParams.alpha * inst.rangeRate + (1-inst.unrollingParams.alpha) * instanteneousRangeRate;
    uVec = gtrack_unrollRadialVelocity(inst.maxRadialVelocity, inst.rangeRate, rvIn);
    
    rvError = (inst.H_s.array.doppler - uVec)/uVec;
    if(abs(rvError) < 0.1)
        inst.velocityHandling = 3;
    end
else
    uVec = gtrack_unrollRadialVelocity(inst.maxRadialVelocity, inst.H_s.array.doppler, uVec);
end
end