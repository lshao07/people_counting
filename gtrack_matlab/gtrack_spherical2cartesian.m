function cart = gtrack_spherical2cartesian(sph)
    range = sph.range;
    azim = sph.azimuth;
    elev = sph.elevation;
    doppler = sph.doppler;
    cart(1) = range*cos(elev)*sin(azim);
    cart(2) = range*cos(elev)*cos(azim);
    cart(3) = range*sin(elev);
    
    cart(4) = doppler*cos(elev)*sin(azim);
    cart(5) = doppler*cos(elev)*cos(azim);
    cart(6) = doppler*sin(elev);
    
    cart(7) = 0;
    cart(8) = 0;
    cart(9) = 0;
end