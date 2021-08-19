function  handle = DPU_TrackerProc_config()
%     global trackerProcObj;
    
    cfgFileName = 'D:\matlab_workspace\people_count\ODS_6m_default.cfg';
    radarCfg = readRadarCfg(cfgFileName);

    %/* Convinience pointer to config params */
%     pStaticCfg = pConfigIn.staticCfg;
    
    %     %/*! @brief      Application level parameters */ 
%     inst.staticCfg.trackerEnabled;         
    %/*! @brief      Application level parameters */ 
    pStaticCfg.staticCfg.sensorAzimuthTilt = radarCfg.sensorPosition(2)*3.14/180;
%     inst.staticCfg.allocationParams;
%     %/*! @brief      Application level parameters */ 
%     inst.staticCfg.gatingParams;
%     %/*! @brief      Application level parameters */ 
%     inst.staticCfg.stateParams;
%     %/*! @brief      Application level parameters */ 
%     inst.staticCfg.variationParams;
    %/*! @brief      Application level parameters */ 
    pStaticCfg.staticCfg.sceneryParams.sensorPosition.x = 0;
    pStaticCfg.staticCfg.sceneryParams.sensorPosition.y = 0;
%     pStaticCfg.staticCfg.sceneryParams.sensorPosition.z = radarCfg.sensorPosition(1);
%     %/*! @brief      Application level parameters */ 
%     inst.staticCfg.accelerationParams[3];
%     %/*! @brief      Application level parameters */
%     inst.staticCfg.trackingParamSet;       
    
    pStaticCfg.staticCfg.gtrackModuleConfig.maxRadialVelocity = radarCfg.trackingCfg(5);%
    pStaticCfg.staticCfg.gtrackModuleConfig.radialVelocityResolution = radarCfg.trackingCfg(6);%
    pStaticCfg.staticCfg.gtrackModuleConfig.initialRadialVelocity = 0; %For people counting %
    pStaticCfg.staticCfg.gtrackModuleConfig.deltaT = radarCfg.trackingCfg(7)*0.001;%
    pStaticCfg.staticCfg.gtrackModuleConfig.maxAcceleration = radarCfg.maxAcceleration;%
    
    defaultSceneryParams.numBoundaryBoxes = 1;
    defaultSceneryParams.boundaryBox = radarCfg.boundaryBox;
    defaultSceneryParams.numStaticBoxes = 1;
    defaultSceneryParams.staticBox = radarCfg.staticBoundaryBox;

    defaultGatingParams.gain = radarCfg.gatingParam(1);
    defaultGatingParams.limits.width = radarCfg.gatingParam(2);
    defaultGatingParams.limits.depth = radarCfg.gatingParam(3);
    defaultGatingParams.limits.height = radarCfg.gatingParam(4);
    defaultGatingParams.limits.vel = radarCfg.gatingParam(5);

    defaultStateParams.det2actThre = radarCfg.stateParam(2);
    defaultStateParams.det2freeThre = radarCfg.stateParam(3);
    defaultStateParams.active2freeThre = radarCfg.stateParam(4);
    defaultStateParams.static2freeThre = radarCfg.stateParam(5);
    defaultStateParams.exit2freeThre = radarCfg.stateParam(6);

    defaultAllocationParams.snrThre = radarCfg.allocationParam(1);
    defaultAllocationParams.snrThreObscured = radarCfg.allocationParam(2);
    defaultAllocationParams.velocityThre = radarCfg.allocationParam(3);
    
    defaultAllocationParams.pointsThre = radarCfg.allocationParam(4);
%     defaultAllocationParams.pointsThre = 10;
    
    defaultAllocationParams.maxDistanceThre = radarCfg.allocationParam(5);
    defaultAllocationParams.maxVelThre = radarCfg.allocationParam(6);

    defaultVariationParams.widthStd = 1/3.46;
    defaultVariationParams.depthStd = 1/3.46;
    defaultVariationParams.heightStd = 1/3.46;
    defaultVariationParams.dopplerStd = 2;

%     gtrackModuleConfig = struct();
%     gtrackModuleConfig.gatingParams = defaultGatingParams;%
%     gtrackModuleConfig.stateParams = defaultStateParams;%
%     gtrackModuleConfig.allocationParams = defaultAllocationParams;%
%     gtrackModuleConfig.variationParams = defaultVariationParams;%
%     gtrackModuleConfig.sceneryParams = defaultSceneryParams;%
    
    pStaticCfg.staticCfg.gatingParams = defaultGatingParams;%
    pStaticCfg.staticCfg.stateParams = defaultStateParams;%
    pStaticCfg.staticCfg.allocationParams = defaultAllocationParams;%
    pStaticCfg.staticCfg.variationParams = defaultVariationParams;%
    pStaticCfg.staticCfg.sceneryParams = defaultSceneryParams;%
    pStaticCfg.staticCfg.accelerationParams = radarCfg.maxAcceleration;%
    
    defaultUnrollingParams.alpha = 0.5;
    defaultUnrollingParams.confidence = 0.1;
    
    pStaticCfg.staticCfg.gtrackModuleConfig.advParams.gatingParams = defaultGatingParams;  
    pStaticCfg.staticCfg.gtrackModuleConfig.advParams.stateParams = defaultStateParams;%
    pStaticCfg.staticCfg.gtrackModuleConfig.advParams.allocationParams = defaultAllocationParams;%
    pStaticCfg.staticCfg.gtrackModuleConfig.advParams.variationParams = defaultVariationParams;%
    pStaticCfg.staticCfg.gtrackModuleConfig.advParams.sceneryParams = defaultSceneryParams;%
    pStaticCfg.staticCfg.gtrackModuleConfig.advParams.unrollingParams = defaultUnrollingParams;
    
    pStaticCfg.staticCfg.gtrackModuleConfig.maxNumPoints = radarCfg.trackingCfg(3);
    pStaticCfg.staticCfg.gtrackModuleConfig.maxNumTracks = radarCfg.trackingCfg(4);
    
    trackerProcObj.pointCloudSize =  pStaticCfg.staticCfg.gtrackModuleConfig.maxNumPoints;
    
    trackerProcObj.targetDescrHandle.currentDescr = 0;
    
    targetList = [];
    for n = 1:pStaticCfg.staticCfg.gtrackModuleConfig.maxNumTracks
        targetList = cat(2, targetList, struct());
    end
    
    trackerProcObj.targetDescrHandle.tList = cat(1, targetList, targetList);

    trackerProcObj.targetDescrHandle.tIndex = cat(1, zeros(1, pStaticCfg.staticCfg.gtrackModuleConfig.maxNumPoints), zeros(1, pStaticCfg.staticCfg.gtrackModuleConfig.maxNumPoints));
    
    trackerProcObj.pDpuCfg.staticCfg = pStaticCfg.staticCfg;
    
    trackerProcObj.gtrackHandle = gtrack_create(pStaticCfg.staticCfg.gtrackModuleConfig);
    
    handle = trackerProcObj;

    %     DPU_trackerProc_updateParamTables(pStaticCfg);

end