function rf = Bot_RF(meshInfo,spInfo,eps)
%%%===========================Copyright==================================%%%
	%%%   Version Nov. 2024
	%%%
	%%%   Gengxuan Zhu <zhugx@zju.edu.cn>
	%%%   Institute of Applied Mechanics,Zhejiang University
	%%%
	%%%===========================Description================================%%%
	%%% This is a function to calculate coefficient matrix for bottom edge.
	%%%======================================================================%%%
    
    % Bot edge
    sptEl = spInfo.sptEl;
    xEle = meshInfo.xEle;
    nelx = meshInfo.nelx;
    nspx = spInfo.nspx;
    % Get bottom edge strain points
    bot_spsEl = sptEl(1:nelx,1:nspx);
    bot_sptEl = -ones(nelx,nspx);    
    % Get fixed element points
    xEle = num2cell(xEle(1:nelx,:),1);
    % Calculate shape function values
    % N1 N2 J22
    N1 = (1 - bot_spsEl) .* (1 - bot_sptEl) ./ 4;
    N2 = (1 + bot_spsEl) .* (1 - bot_sptEl) ./ 4;
    dNds = {-(1 - bot_sptEl) / 4,(1 - bot_sptEl) / 4,(1 + bot_sptEl) / 4,-(1 + bot_sptEl) / 4};
    dxds = cellfun(@(dNds,xEle) dNds .* xEle,dNds,xEle,'UniformOutput',false);
    dxds = dxds{1} + dxds{2} + dxds{3} + dxds{4};
    % Get strain data
    eps22 = eps{4};
    eps22 = eps22(1:nelx,1:nspx);
    eps12 = eps{2};
    eps12 = eps12(1:nelx,1:nspx);
    % Node 1
    rf11_sp = -eps12 .* N1 .* N1 .* dxds;    
    rf12_sp = -eps12 .* N1 .* N2 .* dxds;    
    rf21_sp = -eps22 .* N1 .* N1 .* dxds;
    rf22_sp = -eps22 .* N1 .* N2 .* dxds;

    rf11 = Trapezoidal(bot_sptEl,rf11_sp);
    rf12 = Trapezoidal(bot_sptEl,rf12_sp);
    rf21 = Trapezoidal(bot_sptEl,rf21_sp);
    rf22 = Trapezoidal(bot_sptEl,rf22_sp);
    % Node 4
    rf31_sp = -eps12 .* N1 .* N2 .* dxds;    
    rf32_sp = -eps12 .* N2 .* N2 .* dxds;    
    rf41_sp = -eps22 .* N1 .* N2 .* dxds;    
    rf42_sp = -eps22 .* N2 .* N2 .* dxds;
    
    rf31 = Trapezoidal(bot_sptEl,rf31_sp);
    rf32 = Trapezoidal(bot_sptEl,rf32_sp);
    rf41 = Trapezoidal(bot_sptEl,rf41_sp);
    rf42 = Trapezoidal(bot_sptEl,rf42_sp);
    rf = cat(2,rf11,rf21,rf31,rf41,rf12,rf22,rf32,rf42);
end