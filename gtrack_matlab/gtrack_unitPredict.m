function gtrack_unitPredict(uid)
    global trackerProcObj;
    trackerProcObj.gtrackHandle.hTrack(uid).heartBeatCount = trackerProcObj.gtrackHandle.hTrack(uid).heartBeatCount + 1;
    if ~trackerProcObj.gtrackHandle.hTrack(uid).isTargetStatic
        trackerProcObj.gtrackHandle.hTrack(uid).S_apriori_hat = (trackerProcObj.gtrackHandle.hTrack(uid).F*trackerProcObj.gtrackHandle.hTrack(uid).S_hat').';
        P_apr_tmp = trackerProcObj.gtrackHandle.hTrack(uid).F*trackerProcObj.gtrackHandle.hTrack(uid).P_hat*trackerProcObj.gtrackHandle.hTrack(uid).F'+ trackerProcObj.gtrackHandle.hTrack(uid).Q;
        trackerProcObj.gtrackHandle.hTrack(uid).P_apriori_hat = (P_apr_tmp + P_apr_tmp.')/2;%make matrix symmetric
    else
        trackerProcObj.gtrackHandle.hTrack(uid).S_apriori_hat = trackerProcObj.gtrackHandle.hTrack(uid).S_hat;
        trackerProcObj.gtrackHandle.hTrack(uid).P_apriori_hat = trackerProcObj.gtrackHandle.hTrack(uid).P_hat;
    end
    trackerProcObj.gtrackHandle.hTrack(uid).H_s.array = gtrack_carttosph(trackerProcObj.gtrackHandle.hTrack(uid).S_apriori_hat).';
    x_limit = trackerProcObj.gtrackHandle.hTrack(uid).gatingParams.limits.depth/2;
    el_limit = trackerProcObj.gtrackHandle.hTrack(uid).gatingParams.limits.width/2/trackerProcObj.gtrackHandle.hTrack(uid).H_s.array.range;
    az_limit = trackerProcObj.gtrackHandle.hTrack(uid).gatingParams.limits.height/2/trackerProcObj.gtrackHandle.hTrack(uid).H_s.array.range;
    v_limit = trackerProcObj.gtrackHandle.hTrack(uid).gatingParams.limits.vel/2;
    %update limits
    trackerProcObj.gtrackHandle.hTrack(uid).H_limits.vector.range = x_limit;
    trackerProcObj.gtrackHandle.hTrack(uid).H_limits.vector.elevation = el_limit;
    trackerProcObj.gtrackHandle.hTrack(uid).H_limits.vector.azimuth = az_limit;
    trackerProcObj.gtrackHandle.hTrack(uid).H_limits.vector.doppler = v_limit;
end