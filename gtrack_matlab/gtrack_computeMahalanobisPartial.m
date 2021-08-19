function mdp = gtrack_computeMahalanobisPartial(v,D)
    mdp = v(1:3)*D(1:3,1:3)*v(1:3)';
end