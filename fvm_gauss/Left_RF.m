function rf = Left_RF(meshInfo,spInfo,eps)
%%%===========================Copyright==================================%%%
	%%%   Version Nov. 2024
	%%%
	%%%   Gengxuan Zhu <zhugx@zju.edu.cn>
	%%%   Institute of Applied Mechanics,Zhejiang University
	%%%
	%%%===========================Description================================%%%
	%%% This is a function to calculate coefficient matrix for left edge.
	%%%======================================================================%%%
    
    sptEl = spInfo.sptEl;
    yEle = meshInfo.yEle;
    nelx = meshInfo.nelx;
    nely = meshInfo.nely;
    nspx = spInfo.nspx;
    nspy = spInfo.nspy;
    % Get left edge strain points
    left_spsEl = -ones(nely,nspy);
    left_sptEl = sptEl(1:nelx:end,1:nspx:end);
    % Get fixed element points
    yEle = num2cell(yEle(1:nelx:end,:),1);
    % Calculate shape function values
    % N1 N4 J22
    N1 = (1 - left_spsEl) .* (1 - left_sptEl) ./ 4;
    N4 = (1 - left_spsEl) .* (1 + left_sptEl) ./ 4;
    dNdt = {-(1 - left_spsEl) / 4,-(1 + left_spsEl) / 4,(1 + left_spsEl) / 4,(1 - left_spsEl) / 4};
    dydt = cellfun(@(dNdt,yEle) dNdt .* yEle,dNdt,yEle,'UniformOutput',false);
    dydt = dydt{1} + dydt{2} + dydt{3} + dydt{4};
    % Get strain data
    eps11 = eps{1};
    eps11 = eps11(1:nelx:end,1:nspx:end);
    eps12 = eps{2};
    eps12 = eps12(1:nelx:end,1:nspx:end);
    % Node 1
    rf11_sp = -eps11 .* N1 .* N1 .* dydt;    
    rf12_sp = -eps11 .* N1 .* N4 .* dydt;    
    rf21_sp = -eps12 .* N1 .* N1 .* dydt;    
    rf22_sp = -eps12 .* N1 .* N4 .* dydt;
    
    rf11 = Trapezoidal(left_sptEl,rf11_sp);
    rf12 = Trapezoidal(left_sptEl,rf12_sp);
    rf21 = Trapezoidal(left_sptEl,rf21_sp);
    rf22 = Trapezoidal(left_sptEl,rf22_sp);
    % Node 4
    rf31_sp = -eps11 .* N1 .* N4 .* dydt;    
    rf32_sp = -eps11 .* N4 .* N4 .* dydt;    
    rf41_sp = -eps12 .* N1 .* N4 .* dydt;    
    rf42_sp = -eps12 .* N4 .* N4 .* dydt;
    
    rf31 = Trapezoidal(left_sptEl,rf31_sp);
    rf32 = Trapezoidal(left_sptEl,rf32_sp);
    rf41 = Trapezoidal(left_sptEl,rf41_sp);
    rf42 = Trapezoidal(left_sptEl,rf42_sp);
    rf = cat(2,rf11,rf21,rf31,rf41,rf12,rf22,rf32,rf42);
end