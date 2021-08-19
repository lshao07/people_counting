function limits = gtrack_calcMeasurementLimits(range,gate_limits)
    FLT_MIN = 1.175494351e-38;
    FLT_MAX = 3.402823466e+38;
    if(gate_limits.depth <= FLT_MIN)
		limits.range = FLT_MAX;
	else
		limits.range = gate_limits.depth/2;
    end
    if(gate_limits.height <= FLT_MIN)
		limits.elevation = FLT_MAX;
	else
		limits.elevation = atan((gate_limits.height/2)/range);
    end
	if(gate_limits.width <= FLT_MIN)
		limits.azimuth = FLT_MAX;
	else
		limits.azimuth = atan((gate_limits.width/2)/range);
    end
	if(gate_limits.vel <= FLT_MIN)
		limits.doppler = FLT_MAX;
	else
		limits.doppler = gate_limits.vel/2;
    end
end