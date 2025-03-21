function [eps11,eps22,eps12] = Test_Strain_Data(spInfo)
%%%===========================Copyright==================================%%%
	%%%   Version Oct. 2024
	%%%
	%%%   Gengxuan Zhu <zhugx@zju.edu.cn>
	%%%   Institute of Applied Mechanics,Zhejiang University
	%%%
	%%%===========================Description================================%%%
	%%% This is a function to generate testing strain data.
	%%% ----INPUT
	%%% spInfo.nstrain & spnEl: <scalar> = number of strain points and per element
    %%% spInfo.eleSpID: <vector,[Num of element , nstrain]> = strain point ID for each element
    %%% spInfo.spxEl & spyEl: <vector,[Num of element , nstrain]> = X & Y coordinates for each element 
	%%%
	%%% ----OUTPUT
	%%% eps11: <vector,[Num of Elements , spInfo.spnEl]> = strain data eps11 of each element
	%%% eps22: <vector,[Num of Elements , spInfo.spnEl]> = strain data eps22 of each element
	%%% eps12: <vector,[Num of Elements , spInfo.spnEl]> = strain data eps12 of each element
	%%%======================================================================%%%
    
    E11 = ones(1,spInfo.nstrain);
    E22 = ones(1,spInfo.nstrain);
    E12 = 0.5 * ones(1,spInfo.nstrain);
    eps11 = E11(spInfo.eleSpID);
    eps22 = E22(spInfo.eleSpID);
    eps12 = E12(spInfo.eleSpID);
end