function [spsEl,sptEl,w] = Gauss_Point(nEl)
%%%===========================Copyright==================================%%%
	%%%   Version Nov. 2024
	%%%
	%%%   Gengxuan Zhu <zhugx@zju.edu.cn>
	%%%   Institute of Applied Mechanics,Zhejiang University
	%%%
	%%%===========================Description================================%%%
	%%% This is a function to calculate local coordinates of gauss points.
	%%%======================================================================%%%
    gpt = sqrt(3/5);
    s = [-gpt,0,-gpt,-gpt,0,gpt,-gpt,0,gpt];
	t = [-gpt,-gpt,-gpt,0,0,0,gpt,gpt,gpt];
    w = [25/81,40/81,25/81,40/81,64/81,40/81,25/81,40/81,25/81];
    spsEl = repmat(s,nEl,1);
    sptEl = repmat(t,nEl,1);
end