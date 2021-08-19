function c = gtrack_sph2cart(v)

    [sinAzim, cosAzim] = gtrack_sincosd(v.azimuth);
    [sinElev, cosElev] = gtrack_sincosd(v.elevation);

    c.posX = v.range*cosElev*sinAzim;
    c.posY = v.range*cosElev*cosAzim;
    c.posZ = v.range*sinElev;
end

%{
void gtrack_sph2cart(GTRACK_measurement_vector *v, GTRACK_cartesian_position *c)
{
    float sinAzim, cosAzim;
    float sinElev, cosElev;

    gtrack_sincosd(v->azimuth*RAD2DEG,&sinAzim, &cosAzim);
    gtrack_sincosd(v->elev*RAD2DEG,&sinElev, &cosElev);

    c->posX = v->range*cosElev*sinAzim;
    c->posY = v->range*cosElev*cosAzim;
    c->posZ = v->range*sinElev;
}
%}