function gtrack_unitStart(uid,timeStamp,tid,uCenter)
global trackerProcObj;
VELOCITY_INIT = 0;
TRACK_STATE_DETECTION = 2;

trackerProcObj.gtrackHandle.hTrack(uid).tid = tid;
trackerProcObj.gtrackHandle.hTrack(uid).heartBeatCount = timeStamp;
trackerProcObj.gtrackHandle.hTrack(uid).allocationTime = timeStamp;
trackerProcObj.gtrackHandle.hTrack(uid).allocationRange = uCenter.range;
trackerProcObj.gtrackHandle.hTrack(uid).allocationVelocity = uCenter.doppler;
trackerProcObj.gtrackHandle.hTrack(uid).estNumOfPoints = 0;

% Initialize state and counters 
trackerProcObj.gtrackHandle.hTrack(uid).state = TRACK_STATE_DETECTION;
trackerProcObj.gtrackHandle.hTrack(uid).active2freeCount = 0;
trackerProcObj.gtrackHandle.hTrack(uid).detect2activeCount = 0;
trackerProcObj.gtrackHandle.hTrack(uid).detect2freeCount = 0;
trackerProcObj.gtrackHandle.hTrack(uid).currentStateVectorType = trackerProcObj.gtrackHandle.hTrack(uid).stateVectorType;

trackerProcObj.gtrackHandle.hTrack(uid).isTargetStatic = false;

% Radial Velocity initialization 
% Radial Velocity handling is set to start with range rate filtering 
trackerProcObj.gtrackHandle.hTrack(uid).velocityHandling = VELOCITY_INIT;

u.vector = uCenter;

u.vector.doppler = gtrack_unrollRadialVelocity(trackerProcObj.gtrackHandle.hTrack(uid).maxRadialVelocity, trackerProcObj.gtrackHandle.hTrack(uid).initialRadialVelocity, uCenter.doppler);
trackerProcObj.gtrackHandle.hTrack(uid).rangeRate = u.vector.doppler;

% Initialize a-priori State information */
trackerProcObj.gtrackHandle.hTrack(uid).S_apriori_hat = gtrack_spherical2cartesian(u.vector);
trackerProcObj.gtrackHandle.hTrack(uid).H_s.array = u.vector; % Initialize Hs to measurment vector */

% Initialize measurments error limits */
trackerProcObj.gtrackHandle.hTrack(uid).H_limits.vector = gtrack_calcMeasurementLimits(u.vector.range, trackerProcObj.gtrackHandle.hTrack(uid).gatingParams.limits);

% Initialize a-priori Process covariance */

% P_apriori_hat = eye(6) */
% memcpy(trackerProcObj.gtrackHandle.hTrack(uid).P_apriori_hat, eye6x6, sizeof(eye6x6)); */

% P_apriori_hat = diag([0,0,0.5,0.5,1,1]) */
trackerProcObj.gtrackHandle.hTrack(uid).P_apriori_hat = diag([0,0,0,0.5,0.5,0.5,1,1,1]);
    
trackerProcObj.gtrackHandle.hTrack(uid).gD = zeros(4,4);
trackerProcObj.gtrackHandle.hTrack(uid).estSpread.array = [1.0, 2*pi/180, 10*pi/180, 1.0];
        
trackerProcObj.gtrackHandle.hTrack(uid).G = trackerProcObj.gtrackHandle.hTrack(uid).gatingParams.gain;
trackerProcObj.gtrackHandle.hTrack(uid).sFactor = 1.0;

end