function estDim = gtrack_calcDim(mSpread,R)
    estDim(1) = mSpread(1);
    estDim(2) = 2*R*tan(mSpread(2)/2);
    estDim(3) = 2*R*tan(mSpread(3)/2);
    estDim(4) = mSpread(4);
end