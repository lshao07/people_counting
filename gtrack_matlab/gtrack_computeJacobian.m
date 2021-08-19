function J = gtrack_computeJacobian(S_apr)
x = S_apr(1);
y = S_apr(2);
z = S_apr(3);
v_x = S_apr(4);
v_y = S_apr(5);
v_z = S_apr(6);
r = sqrt(x^2+y^2+z^2);
J = [x/r , y/r ,z/r,0,0,0,0,0,0;
    y/(x^2+y^2) , -x/(x^2+y^2),0,0,0,0,0,0,0;
    -(x/r^2) * (z/sqrt(x^2+y^2)), -(y/r^2) * (z/sqrt(x^2+y^2)), sqrt(x^2+y^2)/r^2,0,0,0,0,0,0;
    (y*(v_x*y - v_y*x) + z*(v_x*z-v_z*x))/r^3,(x*(v_y*x - v_x*y) + z*(v_y*z-v_z*y))/r^3,(x*(v_z*x - v_x*z) + y*(v_z*y-v_y*z))/r^3,x/r,y/r,z/r,0,0,0];
end