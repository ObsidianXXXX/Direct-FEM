function rf = Up_RF(meshInfo,spInfo,eps)
%%%===========================Copyright==================================%%%
	%%%   Version Nov. 2024
	%%%
	%%%   Gengxuan Zhu <zhugx@zju.edu.cn>
	%%%   Institute of Applied Mechanics,Zhejiang University
	%%%
	%%%===========================Description================================%%%
	%%% This is a function to calculate coefficient matrix for up edge.
	%%%======================================================================%%%
    
    % Up edge
    spsEl = spInfo.spsEl;
    sptEl = spInfo.sptEl;
    xEle = meshInfo.xEle;
    nelx = meshInfo.nelx;
    nspx = spInfo.nspx;
    % Get up edge strain points
    up_spsEl = spsEl((end - nelx + 1):end,(end - nspx + 1):end);
    up_sptEl = sptEl((end - nelx + 1):end,(end - nspx + 1):end);
    % Get fixed element points
    xEle = num2cell(xEle((end - nelx + 1):end,:),1);
    % Calculate shape function values
    % N3 N4 J22
    N3 = (1 + up_spsEl) .* (1 + up_sptEl) ./ 4;
    N4 = (1 - up_spsEl) .* (1 + up_sptEl) ./ 4;
    dNds = {-(1 - up_sptEl) / 4,(1 - up_sptEl) / 4,(1 + up_sptEl) / 4,-(1 + up_sptEl) / 4};
    dxds = cellfun(@(dNds,xEle) dNds .* xEle,dNds,xEle,'UniformOutput',false);
    dxds = dxds{1} + dxds{2} + dxds{3} + dxds{4};
    % Get strain data
    eps22 = eps{4};
    eps22 = eps22((end - nelx + 1):end,(end - nspx + 1):end);
    eps12 = eps{2};
    eps12 = eps12((end - nelx + 1):end,(end - nspx + 1):end);
    % Node 3
    rf11_sp = -eps12 .* N3 .* N3 .* dxds;
    rf11 = Trapezoidal(up_sptEl,rf11_sp);
    rf12_sp = -eps12 .* N3 .* N4 .* dxds;
    rf12 = Trapezoidal(up_sptEl,rf12_sp);
    rf21_sp = -eps22 .* N3 .* N3 .* dxds;
    rf21 = Trapezoidal(up_sptEl,rf21_sp);
    rf22_sp = -eps22 .* N3 .* N4 .* dxds;
    rf22 = Trapezoidal(up_sptEl,rf22_sp);
    % Node 4
    rf31_sp = -eps12 .* N3 .* N4 .* dxds;
    rf31 = Trapezoidal(up_sptEl,rf31_sp);
    rf32_sp = -eps12 .* N4 .* N4 .* dxds;
    rf32 = Trapezoidal(up_sptEl,rf32_sp);
    rf41_sp = -eps22 .* N3 .* N4 .* dxds;
    rf41 = Trapezoidal(up_sptEl,rf41_sp);
    rf42_sp = -eps22 .* N4 .* N4 .* dxds;
    rf42 = Trapezoidal(up_sptEl,rf42_sp);
    rf = cat(2,rf11,rf21,rf31,rf41,rf12,rf22,rf32,rf42);
    % Assemble
    %Ivar = meshInfo.upIvar;
    %ndof = meshInfo.ndof;
    %nNod = meshInfo.nNod;
    %rf_C = fsparse(Ivar(:,1),Ivar(:,2),reshape(rf',[],1),[ndof,nNod]);
    %C = C - rf_C;
end