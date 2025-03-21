function Fy = Load_Fy(meshInfo,spInfo,eps,E_load)
%%%===========================Copyright==================================%%%
	%%%   Version Nov. 2024
	%%%
	%%%   Gengxuan Zhu <zhugx@zju.edu.cn>
	%%%   Institute of Applied Mechanics,Zhejiang University
	%%%
	%%%===========================Description================================%%%
	%%% This is a function to calculate resultant force along Y direction.
	%%%======================================================================%%%
    
    spsEl = spInfo.spsEl;
    sptEl = spInfo.sptEl;
    yEle = meshInfo.yEle;
    nelx = meshInfo.nelx;
    nely = meshInfo.nely;
    nspx = spInfo.nspx;
    % Get fixed strain points
    right_spsEl = spsEl(nelx:nelx:end,nspx:nspx:end);
    right_sptEl = sptEl(nelx:nelx:end,nspx:nspx:end);
    % Get fixed element points
    yEle = num2cell(yEle(nelx:nelx:end,:),1);
    % Calculate shape function values
    % N2 N3 J22
    N2 = (1 + right_spsEl) .* (1 - right_sptEl) ./ 4;
    N3 = (1 + right_spsEl) .* (1 + right_sptEl) ./ 4;
    dNdt = {-(1 - right_spsEl) / 4,-(1 + right_spsEl) / 4,(1 + right_spsEl) / 4,(1 - right_spsEl) / 4};
    J22 = cellfun(@(dNdt,yEle) dNdt .* yEle,dNdt,yEle,'UniformOutput',false);
    J22 = J22{1} + J22{2} + J22{3} + J22{4};
    % Get strain data
    eps12 = eps{2};
    eps12 = eps12(nelx:nelx:end,nspx:nspx:end);
    % Calculate force coefficients
    rf1_sp = eps12 .* N2 .* J22;
    rf1 = Trapezoidal(right_sptEl,rf1_sp);
    rf2_sp = eps12 .* N3 .* J22;
    rf2 = Trapezoidal(right_sptEl,rf2_sp);
    rf = cat(2,rf1,rf2);
    % Assemble
    dof_index = [1 1];
    nod_index = [1 2];
	meshnodes = reshape(1:nely,nely,1);
    eleNodesID = meshnodes + [0 1];
    ik = (eleNodesID(:,dof_index))';
    jk = (eleNodesID(:,nod_index))';
	Ivar = [ik(:),jk(:)];
    nody = nely + 1;
    rf_y = fsparse(Ivar(:,1),Ivar(:,2),reshape(rf',[],1),[nody,nody]);
    Fy = sum(rf_y * E_load);
end