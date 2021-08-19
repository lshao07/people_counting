function inst = gtrack_create(config)

inst.maxNumPoints = config.maxNumPoints;
inst.maxNumTracks = config.maxNumTracks;

inst.heartBeat = 0;

defaultUnrollingParams.alpha = 0.5;
defaultUnrollingParams.confidence = 0.1;
params.unrollingParams = defaultUnrollingParams;
params.gatingParams = config.advParams.gatingParams;
params.stateParams = config.advParams.stateParams;
params.allocationParams = config.advParams.allocationParams;
params.variationParams = config.advParams.variationParams;
params.sceneryParams = config.advParams.sceneryParams;

params.deltaT = config.deltaT;
dt = config.deltaT;
dt2 = dt^2;
dt3 = dt^3;
dt4 = dt^4;

params.maxAcceleration = config.maxAcceleration;

% initialize process variance to 1/2 of maximum target acceleration 
varX = (0.5*config.maxAcceleration(1))^2;
varY = (0.5*config.maxAcceleration(2))^2;
varZ = (0.5*config.maxAcceleration(3))^2;

F9 = [			
    1, 0, 0, dt,  0, 0, dt2/2, 0,   0;
    0, 1, 0, 0, dt,  0, 0,   dt2/2, 0;
    0, 0, 1, 0, 0, dt,  0,   0,   dt2/2;
    0, 0, 0, 1, 0, 0, dt,    0,   0;
    0, 0, 0, 0, 1, 0, 0,   dt,    0;
    0, 0, 0, 0, 0, 1, 0,   0,   dt;
    0, 0, 0, 0, 0, 0, 1,   0,   0;
    0, 0, 0, 0, 0, 0, 0,   1,   0;
    0, 0, 0, 0, 0, 0, 0,   0,   1
    ];

Q9 = [ 
    dt4/4*varX,	0,        0,        dt3/2*varX, 0,        0,        dt2/2*varX, 0,        0;
	0,        dt4/4*varY,	0,        0,        dt3/2*varY,	0,        0,        dt2/2*varY, 0;
	0,        0,        dt4/4*varZ,	0,        0,        dt3/2*varZ,	0,        0,        dt2/2*varZ;
	dt3/2*varX,	0,        0,        dt2*varX,	0,        0,        dt*varX,    0,        0;
	0,        dt3/2*varY,	0,        0,        dt2*varY,	0,        0,        dt*varY,    0;
	0,        0,        dt3/2*varZ,	0,        0,        dt2*varZ,	0,        0,        dt*varZ;
	dt2/2*varX,	0,        0,        dt*varX,    0,        0,        1*varX,   0,        0;
	0,        dt2/2*varY,	0,        0,        dt*varY,    0,        0,        1*varY,   0;
    0,        0,        dt2/2*varZ,	0,        0,        dt*varZ,    0,        0,        1*varZ
];

params.F = F9;
params.Q = Q9;

params.maxRadialVelocity = config.maxRadialVelocity;
params.radialVelocityResolution = config.radialVelocityResolution;
params.initialRadialVelocity = config.initialRadialVelocity;


inst.hTrack = []; %hTrack is an array of void pointers
inst.bestScore = [];  %scoreSheet is an array of best scores
inst.bestIndex = []; %association array holds the ids of the best scorers
inst.isUniqueIndex = []; 
inst.allocIndex = []; %allocation array holds the measurement indices of allocation set
inst.uidElem = []; %array of tracking IDs

inst.targetNumTotal = 0;
inst.targetNumCurrent = 0;

inst.freeList.count = 0;
inst.freeList.list = [];

inst.activeList.count = 0;
inst.activeList.list = [];
inst.params = params;

%/* Create unit trackers */
for uid = 1:inst.maxNumTracks
    inst.uidElem(uid) = uid;
    inst.params.uid = uid;
    inst.freeList.list(uid) = inst.uidElem(uid);
    inst.freeList.count = inst.freeList.count + 1;
    inst.hTrack = cat(1, inst.hTrack, gtrack_unitCreate(inst.params));
end




% S = [];
% 
% 
% S = cat(1, S, struct());
% 
% 
% S3DA = [x y z xV yV zV xA yA zA].';
% 
% 
% 
% Jcob3DA = [
%     x/r y/r z/r 0 0 0 0 0 0;
%     y/(x^2+y^2) -x/(x^2+y^2) 0 0 0 0 0 0 0;
%     -(x/r^2)*(z/sqrt(x^2+y^2)) -(y/r^2)*(z/sqrt(x^2+y^2)) sqrt(x^2+y^2)/r^2 0 0 0 0 0 0;
%     y*(xV*y-yV*x)+z*(xV*z-zV*x)/r^3 x*(yV*x-xV*y)+z*(yV*z-zV*y)/r^3 x*(zV*x-xV*z)+y*(zV*y-yV*z)/r^3 x/r y/r z/r 0 0 0
%     ];
% 
% S(1).S3DA = S3DA;
% S(1).F3DA = F3DA;
% S(1).Jcob3DA = Jcob3DA;
end
