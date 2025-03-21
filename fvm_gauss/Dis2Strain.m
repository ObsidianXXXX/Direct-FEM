function [eps11,eps22,eps12] = Dis2Strain(meshInfo,locJ,detJ,dNds,dNdt,ux,uy)
%%%===========================Copyright==================================%%%
	%%%   Version Oct. 2024
	%%%
	%%%   Gengxuan Zhu <zhugx@zju.edu.cn>
	%%%   Institute of Applied Mechanics,Zhejiang University
	%%%
	%%%===========================Description================================%%%
	%%% This is a function to calculate strain values at strain points per
    %%% element from diplacement data.
	%%%======================================================================%%%
    uxEle = num2cell(ux(meshInfo.eleNodsID),1);
    uyEle = num2cell(uy(meshInfo.eleNodsID),1);
    detJ = cell2mat(detJ);
    eps11 = cellfun(@(dNds,dNdt,uxEle) ...
        ((locJ{1} .* dNds + locJ{2} .* dNdt) ./ detJ) .* uxEle, ...
        dNds,dNdt,uxEle,'UniformOutput',false);
    eps22 = cellfun(@(dNds,dNdt,uyEle) ...
        ((locJ{3} .* dNds + locJ{4} .* dNdt) ./ detJ) .* uyEle, ...
        dNds,dNdt,uyEle,'UniformOutput',false);
    eps12 = cellfun(@(dNds,dNdt,uxEle,uyEle) ...
        ((locJ{3} .* dNds + locJ{4} .* dNdt) .* uxEle + ...
         (locJ{1} .* dNds + locJ{2} .* dNdt) .* uyEle) ./ detJ, ...
        dNds,dNdt,uxEle,uyEle,'UniformOutput',false);
    eps11 = eps11{1} + eps11{2} + eps11{3} + eps11{4};
    eps22 = eps22{1} + eps22{2} + eps22{3} + eps22{4};
    eps12 = eps12{1} + eps12{2} + eps12{3} + eps12{4};
end