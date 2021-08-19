function inst = gtrack_unitCreate(params)
    %/* Parameters */
    inst.ec = zeros(4,4);
	inst.gatingParams = params.gatingParams;
	inst.stateParams = params.stateParams;
	inst.allocationParams = params.allocationParams;
	inst.unrollingParams = params.unrollingParams;
	inst.variationParams = params.variationParams;
    inst.sceneryParams = params.sceneryParams;
    
    inst.maxAcceleration = params.maxAcceleration;
    
    inst.uid = params.uid;
    inst.isTargetStatic = false;
	inst.maxRadialVelocity = params.maxRadialVelocity;

	inst.radialVelocityResolution = params.radialVelocityResolution;
% 	inst.verbose = params.verbose;
	inst.initialRadialVelocity = params.initialRadialVelocity;

	inst.F = params.F; 
	inst.Q = params.Q; 
    
    inst.stateVectorType = 'GTRACK_STATE_VECTORS_3DA';
    inst.stateVectorDimNum = 3;
    inst.stateVectorDimLength = 3;
	inst.stateVectorLength = 9;
	inst.measurementVectorLength = 4;
    inst.dt = params.deltaT;
	inst.state = 0;
end