% function [inst,t]  = gtrack_step(handle, points, mNum, mIndex)
function [tNum, targetDescr] = gtrack_step(var,mNum, targetDescr, targetIndex,tNum)
    global trackerProcObj;
    global points;
    trackerProcObj.gtrackHandle.heartBeat = trackerProcObj.gtrackHandle.heartBeat + 1;
    
    if mNum > trackerProcObj.gtrackHandle.maxNumPoints
        mNum = trackerProcObj.gtrackHandle.maxNumPoints;
    end
    
    for n = 1:trackerProcObj.gtrackHandle.maxNumPoints
        trackerProcObj.gtrackHandle.isUniqueIndex(n) = 1;
    end
    
    for n = 1:mNum
        trackerProcObj.gtrackHandle.bestScore(n) = 3.402823466e38; 
        
        if trackerProcObj.gtrackHandle.params.sceneryParams.numBoundaryBoxes
            %/* If boundaries exists, set to outside, and overwrite if inside the boundary */
            trackerProcObj.gtrackHandle.bestIndex(n) = 254;    %GTRACK_ID_POINT_BEHIND_THE_WALL
            
            pos = gtrack_spherical2cartesian(points(n));
            
            for numBoxes = 1:trackerProcObj.gtrackHandle.params.sceneryParams.numBoundaryBoxes
                flag = gtrack_isPointInsideBox(pos, trackerProcObj.gtrackHandle.params.sceneryParams.boundaryBox);
                if  flag
                    %/* Valid points */
                    trackerProcObj.gtrackHandle.bestIndex(n) = 255;    %GTRACK_ID_POINT_NOT_ASSOCIATED
                    break;
                end
            
            end
        else
            %/* Valid points */
		    trackerProcObj.gtrackHandle.bestIndex(n) = 255;    %GTRACK_ID_POINT_NOT_ASSOCIATED
        end
        
    end
    
    %MARK：
    %需要补全module部分代码
	gtrack_modulePredict();
	gtrack_moduleAssociate(mNum);
	gtrack_moduleAllocate(mNum);
	gtrack_moduleUpdate(var, mNum);
	[tNum,targetDescr] = gtrack_moduleReport(tNum,targetDescr);
    
    
    if any(targetIndex)
        for n=1:mNum
            targetIndex(n) = trackerProcObj.gtrackHandle.bestIndex(n);
        end
    end
    %{
    if(bench != NULL)
	{
		/* Collect benchamrks */
		bench[GTRACK_BENCHMARK_SETUP] = gtrack_getCycleCount();
		gtrack_modulePredict(trackerProcObj);
		bench[GTRACK_BENCHMARK_PREDICT] = gtrack_getCycleCount();
		gtrack_moduleAssociate(trackerProcObj, points, mNum);
		bench[GTRACK_BENCHMARK_ASSOCIATE] = gtrack_getCycleCount();
		gtrack_moduleAllocate(trackerProcObj, points, mNum);
		bench[GTRACK_BENCHMARK_ALLOCATE] = gtrack_getCycleCount();
		gtrack_moduleUpdate(trackerProcObj, points, var, mNum);
		bench[GTRACK_BENCHMARK_UPDATE] = gtrack_getCycleCount();
		gtrack_moduleReport(trackerProcObj, t, tNum);
		bench[GTRACK_BENCHMARK_REPORT] = gtrack_getCycleCount();
	}
	else
	{
		gtrack_modulePredict(trackerProcObj);
		gtrack_moduleAssociate(trackerProcObj, points, mNum);
		gtrack_moduleAllocate(trackerProcObj, points, mNum);
		gtrack_moduleUpdate(trackerProcObj, points, var, mNum);
		gtrack_moduleReport(trackerProcObj, t, tNum);
	}

	/* If requested, report uids associated with measurment vector */
	if(mIndex != 0) {
		for(n=0; n< mNum; n++) {
			mIndex[n] = trackerProcObj->bestIndex[n];
		}
	}
    %}
    
    
%     GTRACK_targetDesc.uid;
%     %/**  @brief   Tracking Unit Identifier */
% 	GTRACK_targetDesc.uid;
% 	%/**  @brief   Target Identifier */
% 	GTRACK_targetDesc.tid;
% 	%/**  @brief   State vector */
% 	GTRACK_targetDesc.S[GTRACK_STATE_VECTOR_SIZE];
% 	%/**  @brief   Group covariance matrix */
% 	GTRACK_targetDesc.EC[GTRACK_MEASUREMENT_VECTOR_SIZE*GTRACK_MEASUREMENT_VECTOR_SIZE];
% 	%/**  @brief   Gain factor */
% 	GTRACK_targetDesc.G;
% 	%/**  @brief   Estimated target dimensions: depth, width, [height], doppler */
% 	GTRACK_targetDesc.dim[GTRACK_MEASUREMENT_VECTOR_SIZE];


end