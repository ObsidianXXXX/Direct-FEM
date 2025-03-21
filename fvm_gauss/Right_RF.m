function rf = Right_RF(meshInfo,spInfo,eps)
%%%===========================Copyright==================================%%%
	%%%   Version Nov. 2024
	%%%
	%%%   Gengxuan Zhu <zhugx@zju.edu.cn>
	%%%   Institute of Applied Mechanics,Zhejiang University
	%%%
	%%%===========================Description================================%%%
	%%% This is a function to calculate coefficient matrix for right edge.
	%%%======================================================================%%%
    
    sptEl = spInfo.sptEl;
    yEle = meshInfo.yEle;
    nelx = meshInfo.nelx;
    nely = meshInfo.nely;
    nspx = spInfo.nspx;
    nspy = spInfo.nspy;
    % Get right edge strain points
    right_spsEl = ones(nely,nspy);
    right_sptEl = sptEl(nelx:nelx:end,nspx:nspx:end);
    % Get strain data
    eps11 = eps{1};
    eps11 = eps11(nelx:nelx:end,nspx:nspx:end);
    eps12 = eps{2};
    eps12 = eps12(nelx:nelx:end,nspx:nspx:end);
    % Get right element points
    yEle = num2cell(yEle(nelx:nelx:end,:),1);
    % Calculate shape function values
    % N2 N3 J22
    N2 = (1 + right_spsEl) .* (1 - right_sptEl) ./ 4;
    N3 = (1 + right_spsEl) .* (1 + right_sptEl) ./ 4;
    dNdt = {-(1 - right_spsEl) / 4,-(1 + right_spsEl) / 4,(1 + right_spsEl) / 4,(1 - right_spsEl) / 4};
    dydt = cellfun(@(dNdt,yEle) dNdt .* yEle,dNdt,yEle,'UniformOutput',false);
    dydt = dydt{1} + dydt{2} + dydt{3} + dydt{4};
    
    % Node 2
    rf11_sp = eps11 .* N2 .* N2 .* dydt;
    rf12_sp = eps11 .* N2 .* N3 .* dydt;
    rf21_sp = eps12 .* N2 .* N2 .* dydt;
    rf22_sp = eps12 .* N2 .* N3 .* dydt;
    
    rf11 = Trapezoidal(right_sptEl,rf11_sp);
    rf12 = Trapezoidal(right_sptEl,rf12_sp);
    rf21 = Trapezoidal(right_sptEl,rf21_sp);
    rf22 = Trapezoidal(right_sptEl,rf22_sp);
    
    % Node 3
    rf31_sp = eps11 .* N2 .* N3 .* dydt;
    rf32_sp = eps11 .* N3 .* N3 .* dydt;
    rf41_sp = eps12 .* N2 .* N3 .* dydt;
    rf42_sp = eps12 .* N3 .* N3 .* dydt;
    
    rf31 = Trapezoidal(right_sptEl,rf31_sp);
    rf32 = Trapezoidal(right_sptEl,rf32_sp);
    rf41 = Trapezoidal(right_sptEl,rf41_sp);
    rf42 = Trapezoidal(right_sptEl,rf42_sp);
    
    rf = cat(2,rf11,rf21,rf31,rf41,rf12,rf22,rf32,rf42);
end