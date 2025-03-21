function C = Get_C(meshInfo,spInfo,N,dNds,dNdt,locJ,eps,w)
%%%===========================Copyright==================================%%%
	%%%   Version Nov. 2024
	%%%
	%%%   Gengxuan Zhu <zhugx@zju.edu.cn>
	%%%   Institute of Applied Mechanics,Zhejiang University
	%%%
	%%%===========================Description================================%%%
	%%% This is a function to calculate coefficient matrix for all
    %%% elements.
	%%%======================================================================%%%
    spValue = Sp_Value(meshInfo,spInfo,N,dNds,dNdt,locJ,eps);
    value = cellfun(@(spValue) spValue * w',spValue,'UniformOutput',false);
    vEle = cell2mat(reshape(value,1,[]));
    coe_list = reshape(vEle',[],1);
    Ivar = meshInfo.Ivar;
    ndof = meshInfo.ndof;
    nNod = meshInfo.nNod;
    C = fsparse(Ivar(:,1),Ivar(:,2),coe_list,[ndof,nNod]);
end