function doppler = gtrack_unrollRadialVelocity(maxRadialVelocity,center,current)
    distance = center - current;
    if distance >= 0
        f = floor((distance + maxRadialVelocity)/(2*maxRadialVelocity));
        doppler = current + 2*maxRadialVelocity*f;
    else
        f = floor((maxRadialVelocity-distance)/(2*maxRadialVelocity));
        doppler = current - 2*maxRadialVelocity*f;
    end
end