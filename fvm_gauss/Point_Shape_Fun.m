function [N,dNds,dNdt] = Point_Shape_Fun(zeta)
%%%===========================Copyright==================================%%%
	%%%   Version Oct. 2024
	%%%
	%%%   Gengxuan Zhu <zhugx@zju.edu.cn>
	%%%   Institute of Applied Mechanics,Zhejiang University
	%%%
	%%%===========================Description================================%%%
	%%% This is a function to calculate shape function values at a certain
    %%% point.
	%%%======================================================================%%%
    s = zeta(1);
    t = zeta(2);
    N = [(1-s) * (1-t)/4,(1+s) * (1-t)/4,(1+s) * (1+t)/4,(1-s) * (1+t)/4];
    dNds = [-(1-t)/4,(1-t)/4,(1+t)/4,-(1+t)/4];
    dNdt = [-(1-s)/4,-(1+s)/4,(1+s)/4,(1-s)/4];
end