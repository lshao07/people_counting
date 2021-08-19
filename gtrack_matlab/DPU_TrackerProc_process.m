function targetList = DPU_TrackerProc_process(numObjsIn, detObjIn, detObjInSideInfo)
    
    global trackerProcObj;
    global points;
    global outTrackerProc;

    targetDescr = [];
    tNum = 0;
    
    trackerProcObj.inProgress = true;
    trackerProcObj.pointCloud = [];
    
    pDpuCfg = trackerProcObj.pDpuCfg;
    
    pDpuCfg.res.detObjIn = detObjIn;
    pDpuCfg.res.numDetObjIn = numObjsIn;
    pDpuCfg.res.detObjInSideInfo = detObjInSideInfo;
        
    if (pDpuCfg.res.numDetObjIn) > (pDpuCfg.staticCfg.gtrackModuleConfig.maxNumPoints)
        mNum = pDpuCfg.staticCfg.gtrackModuleConfig.maxNumPoints;
    else
        mNum = pDpuCfg.res.numDetObjIn;
    end
    
    %/* Copy input points to tracker object local point cloud structure */
    for n = 1:mNum
        trackerProcObj.pointCloud(n).range = pDpuCfg.res.detObjIn.range(n);     
        trackerProcObj.pointCloud(n).azimuth = pDpuCfg.res.detObjIn.azim(n) + pDpuCfg.staticCfg.sensorAzimuthTilt;     
        trackerProcObj.pointCloud(n).elevation = pDpuCfg.res.detObjIn.elev(n);     
        trackerProcObj.pointCloud(n).doppler = pDpuCfg.res.detObjIn.doppler(n);     
        trackerProcObj.pointCloud(n).snr = pDpuCfg.res.detObjIn.snr(n);
    end  
    
    points = trackerProcObj.pointCloud;
    
    %Mark：
    %以下三行原处为trackerproc_3d.c：472,生成一个targetList和targetIndex，
    %目前理解是为target输出预留空间
    currentDescr = trackerProcObj.targetDescrHandle.currentDescr;
    targetList = trackerProcObj.targetDescrHandle.tList(currentDescr + 1, :);
    targetIndex = trackerProcObj.targetDescrHandle.tIndex(currentDescr + 1, :);
   
    %Mark：
    %还需要理解gtrack_step函数的输出和输入
    [tNum, targetDescr] = gtrack_step(0, mNum, targetDescr, targetIndex, tNum);
    
    for n=1:tNum
        targetList(n).tid  = targetDescr(n).uid;
        targetList(n).posX = targetDescr(n).S(1);
        targetList(n).posY = targetDescr(n).S(2);
        targetList(n).posZ = targetDescr(n).S(3);
        targetList(n).velX = targetDescr(n).S(4);
        targetList(n).velY = targetDescr(n).S(5);
        targetList(n).velZ = targetDescr(n).S(6);
        targetList(n).accX = targetDescr(n).dim(1);
        targetList(n).accY = targetDescr(n).dim(2);
        targetList(n).accY = targetDescr(n).dim(3);
    end
    
    %/* Fill output parameters */
    if tNum > 0
        outTrackerProc.numTargets = tNum;
    else
        outTrackerProc.numTargets = 0;
    end
    
   %/* Target Indices exist only when we have both points AND targets */
    if (mNum > 0) && (tNum > 0)
        outTrackerProc.numIndices = mNum;
    else
        outTrackerProc.numIndices = 0;
    end
    
    outTrackerProc.currentTargetDesc = currentDescr;
    
    %/* Clear inProgress state */
    trackerProcObj.inProgress = false;
    
end

