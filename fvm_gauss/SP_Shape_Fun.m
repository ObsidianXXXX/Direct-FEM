function [N,dNds,dNdt] = SP_Shape_Fun(spInfo)
%%%===========================Copyright==================================%%%
	%%%   Version Oct. 2024
	%%%
	%%%   Gengxuan Zhu <zhugx@zju.edu.cn>
	%%%   Institute of Applied Mechanics,Zhejiang University
	%%%
	%%%===========================Description================================%%%
	%%% This is a function to calculate shape function values at strain
    %%% points per element.
	%%%======================================================================%%%
    spsEl = spInfo.spsEl;
    sptEl = spInfo.sptEl;
    
    % Shape value
    N1 = (1 - spsEl) .* (1 - sptEl) ./ 4;
    N2 = (1 + spsEl) .* (1 - sptEl) ./ 4;
    N3 = (1 + spsEl) .* (1 + sptEl) ./ 4;
    N4 = (1 - spsEl) .* (1 + sptEl) ./ 4;
    N = {N1,N2,N3,N4};
    % Derivate value
    dNds1 = -(1 - sptEl) ./ 4;
    dNds2 = (1 - sptEl) ./ 4;
    dNds3 = (1 + sptEl) ./ 4;
    dNds4 = -(1 + sptEl) ./ 4;
    dNds = {dNds1,dNds2,dNds3,dNds4};
    
    dNdt1 = -(1 - spsEl) ./ 4;
    dNdt2 = -(1 + spsEl) ./ 4;
    dNdt3 = (1 + spsEl) ./ 4;
    dNdt4 = (1 - spsEl) ./ 4;
    dNdt = {dNdt1,dNdt2,dNdt3,dNdt4};
end