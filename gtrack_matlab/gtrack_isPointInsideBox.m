function val = gtrack_isPointInsideBox(c, box)
    if  (c(1) > box(1)) && (c(1) < box(2)) && (c(2) > box(3)) && (c(2) < box(4)) && (c(3) > box(5)) && (c(3) < box(6))
        val = true;
    else
        val = false;
    end
end

%{
uint8_t gtrack_isPointInsideBox(GTRACK_cartesian_position *c, GTRACK_boundaryBox *box)
{
    if( (c.posX > box.x1) && (c.posX < box.x2) && 
        (c.posY > box.y1) && (c.posY < box.y2) && 
        (c.posZ > box.z1) && (c.posZ < box.z2) )
        return 1U;
    else
        return 0;
}
%}