function state = gtrack_unitUpdate(uid,var,num)
global trackerProcObj;
global points;
GTRACK_SPREAD_ALPHA = 0.01;
GTRACK_MIN_POINTS_TO_UPDATE_DISPERSION = 3;
FLT_MAX = 3.402823466e+38;
MSIZE = 4;
myPointNum = 0;
myUniquePointNum = 0;
centroidSnr = 0;
J = zeros(4,9);
D = zeros(4,4);
u_sum.array = zeros(1,4);
uvar_sum.array = zeros(1,4);
u_max = ones(1,4)*(-FLT_MAX);
u_min = ones(1,4)*FLT_MAX;
targetCanStop = true;
spreadMin = [0.25, pi/180, 2*pi/180, 0.25];
% trackerProcObj.gtrackHandle.hTrack(uid) = trackerProcObj.gtrackHandle.hTrack(uid);
for n = 1:num
    if trackerProcObj.gtrackHandle.bestIndex(n) == trackerProcObj.gtrackHandle.hTrack(uid).uid
        u.vector = points(n);
        % Compute maximums and minimums for each measurement */
        
        if trackerProcObj.gtrackHandle.isUniqueIndex(n) == 1
            fields = fieldnames(u.vector);
            for m = 1: MSIZE
                k = fields(m);
                key = k{1};
                if u.vector.(key) > u_max(m)
                    u_max(m) = u.vector.(key);
                end
                if u.vector.(key) < u_min(m)
                    u_min(m) = u.vector.(key);
                end
            end
            myUniquePointNum = myUniquePointNum + 1;
        end
        
        if any(var)
            uvar.array = var(n);
            uvar_sum.array = uvar_sum.array + uvar.array;
        end
        
        if myPointNum == 0
            rvPilot = u.vector.doppler;
        else
            u.vector.doppler = gtrack_unrollRadialVelocity(trackerProcObj.gtrackHandle.hTrack(uid).maxRadialVelocity, rvPilot, u.vector.doppler);
        end
        
        centroidSnr = centroidSnr + points(n).snr;
        u_sum.array(1) = u_sum.array(1) + u.vector.range*points(n).snr;
        u_sum.array(2) = u_sum.array(2) + u.vector.elevation*points(n).snr;
        u_sum.array(3) = u_sum.array(3) + u.vector.azimuth*points(n).snr;
        u_sum.array(4) = u_sum.array(4) + u.vector.doppler*points(n).snr;
        myPointNum = myPointNum + 1;
    end
end

if myPointNum
    % Update estimated number of points */
    if myPointNum > trackerProcObj.gtrackHandle.hTrack(uid).estNumOfPoints
        trackerProcObj.gtrackHandle.hTrack(uid).estNumOfPoints = myPointNum;
    else
        trackerProcObj.gtrackHandle.hTrack(uid).estNumOfPoints = 0.99 *trackerProcObj.gtrackHandle.hTrack(uid).estNumOfPoints + 0.01 * myPointNum;
    end
    % lower bound for the estimated number of points is allocation threshold */
    if trackerProcObj.gtrackHandle.hTrack(uid).estNumOfPoints < trackerProcObj.gtrackHandle.hTrack(uid).allocationParams.pointsThre
        trackerProcObj.gtrackHandle.hTrack(uid).estNumOfPoints = trackerProcObj.gtrackHandle.hTrack(uid).allocationParams.pointsThre;
    end
    % Compute measurement centroid as SNR-weighted mean of associated points */
    u_centroid.vector = u_sum.array .* 1/centroidSnr;
    
    % Unroll centroid radial velocity based on target state */
    rvIn = u_centroid.vector(4);
    u_centroid.vector(4) = gtrack_velocityStateHandling(uid,u_centroid.vector);
    
    if trackerProcObj.gtrackHandle.hTrack(uid).isTargetStatic
        % Returning from static condition
        % Restore state information
        trackerProcObj.gtrackHandle.hTrack(uid).S_apriori_hat = trackerProcObj.gtrackHandle.hTrack(uid).S_apriori_saved;
        trackerProcObj.gtrackHandle.hTrack(uid).estSpread.vector.doppler = trackerProcObj.gtrackHandle.hTrack(uid).estDopplerSpread_saved;
        trackerProcObj.gtrackHandle.hTrack(uid).isTargetStatic = false;
    end
    %Update scoring powers
    if trackerProcObj.gtrackHandle.hTrack(uid).sFactor < 1
        trackerProcObj.gtrackHandle.hTrack(uid).sFactor = trackerProcObj.gtrackHandle.hTrack(uid).sFactor + 0.01;
    end
    
    %Compute mean measurment variance, if availbale 
    if any(var)
        uvar_mean.array = uvar_sum.array ./myPointNum;
    end
    % Update measurement spread if we have 2 or more unique points */
    if(myUniquePointNum > 1) 
        fields = fieldnames(trackerProcObj.gtrackHandle.hTrack(uid).H_limits.vector);
        for m = 1: MSIZE
            k = fields(m);
            key = k{1};
            
            spread = u_max(m) - u_min(m);
            % Unbiased spread estimation */
            spread = spread*(myUniquePointNum+1)/(myUniquePointNum-1);
            % if spread can't be smaller than measurment error */
            if(spread < spreadMin(m))
                spread = spreadMin(m);
            end
            % if spread can't be larger than configured limits */
            if(spread > 2*trackerProcObj.gtrackHandle.hTrack(uid).H_limits.vector.(key))
                spread = 2*trackerProcObj.gtrackHandle.hTrack(uid).H_limits.vector.(key);
            end
            if(spread > trackerProcObj.gtrackHandle.hTrack(uid).estSpread.array(m))
                trackerProcObj.gtrackHandle.hTrack(uid).estSpread.array(m) = spread;
            else
                trackerProcObj.gtrackHandle.hTrack(uid).estSpread.array(m) = (1.0-GTRACK_SPREAD_ALPHA)*trackerProcObj.gtrackHandle.hTrack(uid).estSpread.array(m) + GTRACK_SPREAD_ALPHA*spread;
            end
        end
        trackerProcObj.gtrackHandle.hTrack(uid).estDim = gtrack_calcDim(trackerProcObj.gtrackHandle.hTrack(uid).estSpread.array, u_centroid.vector(1));
    end
else
    if(trackerProcObj.gtrackHandle.hTrack(uid).isTargetStatic == false) 
        for n = 1:3
            vel = abs(trackerProcObj.gtrackHandle.hTrack(uid).S_hat(n+3));
            if vel > 1
                targetCanStop = false;
                break;
            end
        end
%         for n=1:trackerProcObj.gtrackHandle.hTrack(uid).stateVectorDimNum
%             vel = abs(trackerProcObj.gtrackHandle.hTrack(uid).S_hat(GTRACK_VEL_DIMENSION*trackerProcObj.gtrackHandle.hTrack(uid).stateVectorDimNum + n));
%             if(vel > GTRACK_MIN_STATIC_VELOCITY)
%                 targetCanStop = false;
%                 break;
%             end
%         end
        if(targetCanStop)
            % Target is static
            trackerProcObj.gtrackHandle.hTrack(uid).isTargetStatic = true;

            % Save state information
            trackerProcObj.gtrackHandle.hTrack(uid).S_apriori_saved = trackerProcObj.gtrackHandle.hTrack(uid).S_apriori_hat;

            % Force zero velocity/zero acceleration
            trackerProcObj.gtrackHandle.hTrack(uid).S_apriori_hat(4:9) = 0;

            trackerProcObj.gtrackHandle.hTrack(uid).estDopplerSpread_saved = trackerProcObj.gtrackHandle.hTrack(uid).estSpread.array(4);
            trackerProcObj.gtrackHandle.hTrack(uid).estSpread.array(4) =  trackerProcObj.gtrackHandle.hTrack(uid).estSpread.array(4) + 2*abs(trackerProcObj.gtrackHandle.hTrack(uid).H_s.array.doppler);

            % Keep process covariance equal to predicted
            trackerProcObj.gtrackHandle.hTrack(uid).P_hat = trackerProcObj.gtrackHandle.hTrack(uid).P_apriori_hat;
            trackerProcObj.gtrackHandle.hTrack(uid).sFactor = 0;
        else
            trackerProcObj.gtrackHandle.hTrack(uid).S_apriori_hat(7:9) = 0;
        end
    end
end


if all(var == 0)
    uvar_mean.array =(trackerProcObj.gtrackHandle.hTrack(uid).estSpread.array(:) ./2).^2;
end
Rm = diag(uvar_mean.array);
tmp_points_list = [];
if myPointNum > GTRACK_MIN_POINTS_TO_UPDATE_DISPERSION
    for n = 1:num
        if trackerProcObj.gtrackHandle.bestIndex(n) == trackerProcObj.gtrackHandle.hTrack(uid).uid
            tmp_points_list = cat(1, tmp_points_list, [points(n).range, points(n).azimuth, points(n).elevation, points(n).doppler]);
        end
    end
    D = cov(tmp_points_list);
            
    % Normalize it
    D = D./myPointNum;
    
    % Update persistant group dispersion based on instantaneous D
    % The factor alpha goes from maximum (1.f) at the first allocation down to minimum of 0.1f once the target is observed for the long time */
    alpha = (myPointNum)/trackerProcObj.gtrackHandle.hTrack(uid).estNumOfPoints;
    
    %  trackerProcObj.gtrackHandle.hTrack(uid)->gD = (1-alpha)*trackerProcObj.gtrackHandle.hTrack(uid)->gD + alpha*D
    trackerProcObj.gtrackHandle.hTrack(uid).gD = (1-alpha).*trackerProcObj.gtrackHandle.hTrack(uid).gD + alpha .* D;
end

J = gtrack_computeJacobian(trackerProcObj.gtrackHandle.hTrack(uid).S_apriori_hat);
JPJ = J*trackerProcObj.gtrackHandle.hTrack(uid).P_apriori_hat*J';

if myPointNum
    % Compute centroid measurement noise covariance matrix Rc used for Kalman updates 
    % First term represents the error in measuring the centroid and decreased with the number of measurments
    % Second term represents the centroid unsertanty due to the fact that not all the memebers observed
    % See page 327, 11A.2 of S.Blackman
    % Rc = Rm/num + alpha*unit.gD*eye(mSize);
    alpha = ((trackerProcObj.gtrackHandle.hTrack(uid).estNumOfPoints-myPointNum))/((trackerProcObj.gtrackHandle.hTrack(uid).estNumOfPoints-1)*myPointNum);
    
    Rc = Rm./myPointNum + alpha.*diag(diag(trackerProcObj.gtrackHandle.hTrack(uid).gD))*eye(4);
    
    u_tilda = u_centroid.vector -  [trackerProcObj.gtrackHandle.hTrack(uid).H_s.array.range, trackerProcObj.gtrackHandle.hTrack(uid).H_s.array.elevation, trackerProcObj.gtrackHandle.hTrack(uid).H_s.array.azimuth, trackerProcObj.gtrackHandle.hTrack(uid).H_s.array.doppler];
    
    %Compute centroid covariance cC = [3x6]x[6x6]x[6x3]+[3x3] 
    % cC = J(:,1:mSize) * obj.P_apriori(1:mSize,1:mSize) * J(:,1:mSize)' + Rc 
    cC = JPJ + Rc;
    cC_inv = inv(cC);
    
    %Compute Kalman Gain K[6x3] = P[6x6]xJ[3x6]'xIC_inv[3x3]=[9x4] 
    %K = obj.P_apriori(1:mSize,1:mSize) * J(:,1:mSize)' * iC_inv 
    K = trackerProcObj.gtrackHandle.hTrack(uid).P_apriori_hat*J'*cC_inv;
    
    % State estimation */
	% obj.S_hat(1:mSize) = obj.S_apriori_hat(1:mSize) + K * u_tilda */
    trackerProcObj.gtrackHandle.hTrack(uid).S_hat = trackerProcObj.gtrackHandle.hTrack(uid).S_apriori_hat + (K*u_tilda.').';
            
    % Covariance estimation
    trackerProcObj.gtrackHandle.hTrack(uid).P_hat = trackerProcObj.gtrackHandle.hTrack(uid).P_apriori_hat - K * J * trackerProcObj.gtrackHandle.hTrack(uid).P_apriori_hat;
    
else    
    % Handling of erasures 
    trackerProcObj.gtrackHandle.hTrack(uid).S_hat =  trackerProcObj.gtrackHandle.hTrack(uid).S_apriori_hat;	
    trackerProcObj.gtrackHandle.hTrack(uid).P_hat =  trackerProcObj.gtrackHandle.hTrack(uid).P_apriori_hat;
end

% Compute groupCovariance gC (will be used in gating) 
% We will use ellipsoidal gating, that acounts for the dispersion of the group, target maneuver, and measurement noise 
% gC = gD + JPJ + Rm 

trackerProcObj.gtrackHandle.hTrack(uid).gC = trackerProcObj.gtrackHandle.hTrack(uid).gD + JPJ + Rm;
trackerProcObj.gtrackHandle.hTrack(uid).gC_inv = inv(trackerProcObj.gtrackHandle.hTrack(uid).gC);

gtrack_unitEvent(uid, myPointNum);
state = trackerProcObj.gtrackHandle.hTrack(uid).state;
% trackerProcObj.gtrackHandle.hTrack(uid) = trackerProcObj.gtrackHandle.hTrack(uid);

end