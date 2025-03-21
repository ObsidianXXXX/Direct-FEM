function F = Load_F(meshInfo,N,loadPoint)
%%%===========================Copyright==================================%%%
	%%%   Version Nov. 2024
	%%%
	%%%   Gengxuan Zhu <zhugx@zju.edu.cn>
	%%%   Institute of Applied Mechanics,Zhejiang University
	%%%
	%%%===========================Description================================%%%
	%%% This is a function to calculate point load.
	%%%======================================================================%%%
    alldof = 2*meshInfo.nNod;
    F = sparse(alldof,1);
    loadDof = loadPoint*2-1;
    FMag = N*meshInfo.eh;
    F(loadDof,1) = FMag;
    F(loadDof(1),1) = F(loadDof(1),1) - FMag/2;
    F(loadDof(end),1) = F(loadDof(end),1) - FMag/2;
end