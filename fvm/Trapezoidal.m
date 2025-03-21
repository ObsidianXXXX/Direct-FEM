function value = Trapezoidal(zeta,spValue)
%%%===========================Copyright==================================%%%
	%%%   Version Nov. 2024
	%%%
	%%%   Gengxuan Zhu <zhugx@zju.edu.cn>
	%%%   Institute of Applied Mechanics,Zhejiang University
	%%%
	%%%===========================Description================================%%%
	%%% This is a function to calculate 1D integral by trapezoidal rule.
	%%%======================================================================%%%
    num = size(spValue,2);
    W(:,1) = (zeta(:,2) - zeta(:,1)) / 2;
    W(:,num) = (zeta(:,end) - zeta(:,end - 1)) / 2;
    W(:,2:num - 1) = (zeta(:,3:end) - zeta(:,1:end - 2)) / 2;
    value = sum(W .* spValue,2);
end