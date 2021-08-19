function flag = gtrack_isPointBehindTarget(p,target)
    if abs(p.azimuth - target.azimuth) < 2*pi/180 && abs(p.elevation - target.elevation) < 2*pi/180 && p.range > target.range       
        flag = true;
    else
        flag = false;
    end
end