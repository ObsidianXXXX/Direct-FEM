function eps = Structured_Strain(eps11,eps22,eps12,nu)
%%%===========================Copyright==================================%%%
	%%%   Version Oct. 2024
	%%%
	%%%   Gengxuan Zhu <zhugx@zju.edu.cn>
	%%%   Institute of Applied Mechanics,Zhejiang University
	%%%
	%%%===========================Description================================%%%
	%%% This is a function to generate structured strain data after multiplying by D.
	%%% ----INPUT
	%%% eps11: <vector,[Num of Elements , spInfo.spnEl]> = strain data eps11 of each element
	%%% eps22: <vector,[Num of Elements , spInfo.spnEl]> = strain data eps22 of each element
	%%% eps12: <vector,[Num of Elements , spInfo.spnEl]> = strain data eps12 of each element
	%%% nu: <scalar> = poisson's ratio
    %%%
    %%% ----OUTPUT
	%%% eps: <cell,4x[Num of Elements , spInfo.spnEl]> = strain data of each element
	%%%======================================================================%%%
    eps_11 = (eps11 - nu .* eps22) ./ (1 - nu ^ 2);
    eps_22 = (eps22 - nu .* eps11) ./ (1 - nu ^ 2);
    eps_12 = eps12 ./ (1 + nu) ./ 2;
    eps = {eps_11,eps_12,eps_12,eps_22};
end