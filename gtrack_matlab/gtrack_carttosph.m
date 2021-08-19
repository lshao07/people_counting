function coord_s = gtrack_carttosph(S_apr)
    x = S_apr(1);
    y = S_apr(2);
    z = S_apr(3);
    v_x = S_apr(4);
    v_y = S_apr(5);
    v_z = S_apr(6);
    r = sqrt(x^2+y^2+z^2);
    v_rad = (x*v_x + y*v_y +z*v_z)/r;
    if y > 0
        phi = atan(x/y);
    elseif y == 0
        phi = pi/2;
    else
        phi = atan(x/y)+pi;
    end
    %range, elevation ,azimuth,doppler
    coord_s.range = r;
    coord_s.elevation = atan(z/sqrt(x^2+y^2));
    coord_s.azimuth = phi;
    coord_s.doppler = v_rad;
    
end