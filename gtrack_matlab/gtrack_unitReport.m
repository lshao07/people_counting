function target = gtrack_unitReport(uid,target)
    global trackerProcObj;
    inst = trackerProcObj.gtrackHandle.hTrack(uid);
	target.uid = inst.uid;
	target.tid = inst.tid;
    target.S = inst.S_hat;
	target.EC = inst.ec;
	target.G = inst.G;
    target.dim = inst.estDim;

end