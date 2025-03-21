function [left,right,up,bot] = Edge_RF(meshInfo,spInfo,eps)
%%%===========================Copyright==================================%%%
	%%%   Version Nov. 2024
	%%%
	%%%   Gengxuan Zhu <zhugx@zju.edu.cn>
	%%%   Institute of Applied Mechanics,Zhejiang University
	%%%
	%%%===========================Description================================%%%
	%%% This is a function to calculate reaction force of nodes on edges.
	%%%======================================================================%%%
    left = Left_RF(meshInfo,spInfo,eps);
    up = Up_RF(meshInfo,spInfo,eps);
    bot = Bot_RF(meshInfo,spInfo,eps);
    right = Right_RF(meshInfo,spInfo,eps);
end