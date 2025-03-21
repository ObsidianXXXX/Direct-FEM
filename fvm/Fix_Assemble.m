function C = Fix_Assemble(meshInfo,fix_rf,C)
%%%===========================Copyright==================================%%%
	%%%   Version Nov. 2024
	%%%
	%%%   Gengxuan Zhu <zhugx@zju.edu.cn>
	%%%   Institute of Applied Mechanics,Zhejiang University
	%%%
	%%%===========================Description================================%%%
	%%% This is a function to modify coefficient matrix.
	%%%======================================================================%%%
    % Assemble
    Ivar = meshInfo.leftIvar;
    ndof = meshInfo.ndof;
    nNod = meshInfo.nNod;
    rf_C = fsparse(Ivar(:,1),Ivar(:,2),reshape(fix_rf',[],1),[ndof,nNod]);
    C = C - rf_C;
end