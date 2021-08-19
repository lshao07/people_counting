function gtrack_unitScore(uid,num)
global trackerProcObj;
global points;
GTRACK_ID_POINT_BEHIND_THE_WALL = 254;
VELOCITY_LOCKED = 3;
GTRACK_MEASUREMENT_VECTOR_SIZE = 4;
GTRACK_NUM_TRACKS_MAX = 250;
% trackerProcObj.gtrackHandle.hTrack(uid) = trackerProcObj.gtrackHandle.hTrack(uid);

logdet = log(det(trackerProcObj.gtrackHandle.hTrack(uid).gC));

for n = 1:num
    if trackerProcObj.gtrackHandle.bestIndex(n) == GTRACK_ID_POINT_BEHIND_THE_WALL
        continue;
    end
    %===========================================================
    u_tilda.vector.range = points(n).range - trackerProcObj.gtrackHandle.hTrack(uid).H_s.array.range;
    u_tilda.vector.elevation = points(n).elevation - trackerProcObj.gtrackHandle.hTrack(uid).H_s.array.elevation;
    u_tilda.vector.azimuth = points(n).azimuth - trackerProcObj.gtrackHandle.hTrack(uid).H_s.array.azimuth;
    u_tilda.vector.doppler = points(n).doppler - trackerProcObj.gtrackHandle.hTrack(uid).H_s.array.doppler;
    
    if trackerProcObj.gtrackHandle.hTrack(uid).velocityHandling < VELOCITY_LOCKED
        % Radial velocity estimation is not yet known, unroll based on velocity measured at allocation time
        rvOut = gtrack_unrollRadialVelocity(trackerProcObj.gtrackHandle.hTrack(uid).maxRadialVelocity, trackerProcObj.gtrackHandle.hTrack(uid).allocationVelocity, points(n).doppler);
        u_tilda.vector.doppler = rvOut - trackerProcObj.gtrackHandle.hTrack(uid).allocationVelocity;
    else
        % Radial velocity estimation is known
        rvOut = gtrack_unrollRadialVelocity(trackerProcObj.gtrackHandle.hTrack(uid).maxRadialVelocity, trackerProcObj.gtrackHandle.hTrack(uid).H_s(4), points(n).doppler);
        u_tilda.vector.doppler = rvOut - trackerProcObj.gtrackHandle.hTrack(uid).H_s(4);
    end
    
    % Any points outside the limits is outside the gate */
	isWithinLimits = true;
    fields = fieldnames(trackerProcObj.gtrackHandle.hTrack(uid).H_limits.vector);
    for m=1:GTRACK_MEASUREMENT_VECTOR_SIZE
        k = fields(m);
        key = k{1};
        if abs(u_tilda.vector.(key)) > trackerProcObj.gtrackHandle.hTrack(uid).H_limits.vector.(key)
            isWithinLimits = false;
            break;
        end
    end
    if (isWithinLimits == false)
        continue;
    end
    u_tilda.array = [u_tilda.vector.range,u_tilda.vector.elevation,u_tilda.vector.azimuth,u_tilda.vector.doppler];
    mdp = gtrack_computeMahalanobisPartial(u_tilda.array, trackerProcObj.gtrackHandle.hTrack(uid).gC_inv);
    
    if mdp < trackerProcObj.gtrackHandle.hTrack(uid).G
        if trackerProcObj.gtrackHandle.bestIndex(n) < GTRACK_NUM_TRACKS_MAX
            % This points is not unique. Clear the indication */
            trackerProcObj.gtrackHandle.isUniqueIndex(n) = 0;
        end
        %Scoring 
        md = gtrack_computeMahalanobis(u_tilda.array, trackerProcObj.gtrackHandle.hTrack(uid).gC_inv);
		score = logdet + md;
		if score < trackerProcObj.gtrackHandle.bestScore(n)
            % If we are the best, register our score, and the index 
            trackerProcObj.gtrackHandle.bestScore(n) = score;
			trackerProcObj.gtrackHandle.bestIndex(n) = trackerProcObj.gtrackHandle.hTrack(uid).uid;
			points(n).vector.doppler = rvOut;
        end
    end
end
trackerProcObj.gtrackHandle.hTrack(uid).ec = trackerProcObj.gtrackHandle.hTrack(uid).gC_inv;
% trackerProcObj.gtrackHandle.hTrack(uid) = trackerProcObj.gtrackHandle.hTrack(uid);
end