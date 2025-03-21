function [meshInfo,spInfo] = Mesh_Sp_Info(l,h,nelx,nely,nspx,nspy)
	%%%===========================Copyright==================================%%%
	%%%   Version Oct. 2024
	%%%
	%%%   Gengxuan Zhu <zhugx@zju.edu.cn>
	%%%   Institute of Applied Mechanics,Zhejiang University
	%%%
	%%%===========================Description================================%%%
	%%% This is a function to mesh a rectangular region with uniform element
    %%% and initialize strain information.
	%%% ----INPUT
	%%% l: <scalar> = length of X
	%%% h: <scalar> = height of Y
	%%% nelx: <scalar> = number of elements along X
	%%% nely: <scalar> = number of elements along Y
    %%% nspx: <scalar> = Number of strain data points per element along X
    %%% nspy: <scalar> = Number of strain data points per element along Y
	%%%
	%%% ----OUTPUT
	%%% meshInfo.coor: <vector,[Num of Nodes , 2]> = coordinates of each node
	%%% meshInfo.eleNodsID: <vector,[Num of element , 4]> = node ID for each element
	%%% meshInfo.el & eh: <scalar> = element length & height
	%%% meshInfo.nEl & nNod: <scalar> = number of elements & nods
	%%% meshInfo.xEl & yEl: <vector,[Num of element , 4]> = X & Y coordinates for each element
    %%% spInfo.nstrain & spnEl: <scalar> = number of strain points and per element
    %%% spInfo.eleSpID: <vector,[Num of element , spnEl]> = strain point ID for each element
    %%% spInfo.spxEl & spyEl: <vector,[Num of element , spnEl]> = X & Y coordinates for each element 
	%%% spInfo.plot_x & plot_y: <vector,[Number of strain points , 1]> = X & Y coordinates
    %%%======================================================================%%%
    el = l/nelx;
	eh = h/nely;

	nodsx = nelx+1;
	nodsy = nely+1;
    ndof = 2 * nodsx * nodsy;

	x = linspace(0,l,nodsx);
	y = linspace(0,h,nodsy);

	[plot_x,plot_y] = meshgrid(x,y);
	coor(:,1) = reshape(plot_x',[],1);
	coor(:,2) = reshape(plot_y',[],1);
    X = coor(:,1); Y = coor(:,2);

	%Assemble relationship
	nEl = nelx*nely;
	meshNods = int32(reshape((1:nodsx*nodsy),nodsx,nodsy));
	% 每个单元左下角节点的第一个自由度
	cVec = reshape(2*meshNods(1:end-1,1:end-1)+1,nEl,1);
	% 定义每个单元四个节点的两个自由度并组装
	cMat = cVec+int32([-2,-1,0,1,2*nelx+[2,3,0,1]]);
	eleNodsID = reshape(meshNods(1:end-1,1:end-1),nEl,1)+int32([0,1,nelx+[2,1]]);
	dof_index = repmat([1,2,3,4,5,6,7,8],1,4);
    nod_index = reshape(repmat([1,2,3,4],8,1),1,[]);
	ik = (cMat(:,dof_index))'; jk = (eleNodsID(:,nod_index))';
	Ivar = [ik(:),jk(:)];
    
    % Left Ivar
    leftdof_index = repmat([1,2,7,8],1,2);
    leftnod_index = reshape(repmat([1,4],4,1),1,[]);
	dof_left = (cMat(1:nelx:end,leftdof_index))';
    nod_left = (eleNodsID(1:nelx:end,leftnod_index))';
	leftIvar = [dof_left(:),nod_left(:)];
    
    % Right Ivar
    rightdof_index = repmat([3,4,5,6],1,2);
    rightnod_index = reshape(repmat([2,3],4,1),1,[]);
	dof_right = (cMat(nelx:nelx:end,rightdof_index))';
    nod_right = (eleNodsID(nelx:nelx:end,rightnod_index))';
	rightIvar = [dof_right(:),nod_right(:)];
    
    % Up Ivar
    updof_index = repmat([5,6,7,8],1,2);
    upnod_index = reshape(repmat([3,4],4,1),1,[]);
	dof_up = (cMat((end - nelx + 1):end,updof_index))';
    nod_up = (eleNodsID((end - nelx + 1):end,upnod_index))';
	upIvar = [dof_up(:),nod_up(:)];

    % Bottom Ivar
    botdof_index = repmat([1,2,3,4],1,2);
    botnod_index = reshape(repmat([1,2],4,1),1,[]);
	dof_bot = (cMat(1:nelx,botdof_index))';
    nod_bot = (eleNodsID(1:nelx,botnod_index))';
	botIvar = [dof_bot(:),nod_bot(:)];

	meshInfo.coord = coor; meshInfo.eleNodsID = eleNodsID;
	meshInfo.eleSize=[el,eh];meshInfo.eleDofID = cMat; meshInfo.Ivar = Ivar;
	meshInfo.nEl = nEl; meshInfo.el = el;  meshInfo.eh = eh; meshInfo.nNod = length(coor(:,1));
    meshInfo.xEle = X(eleNodsID); meshInfo.yEle = Y(eleNodsID);
    meshInfo.nelx = nelx; meshInfo.nely = nely;
    meshInfo.ndof = ndof;
    meshInfo.leftIvar = leftIvar;
    meshInfo.rightIvar = rightIvar;
    meshInfo.upIvar = upIvar;
    meshInfo.botIvar = botIvar;
    
    spx = (nspx - 1) * nelx + 1; % number of strain points along X
    spy = (nspy - 1) * nely + 1; % number of strain points along Y
    nstrain = spx * spy; % number of all strain points
    spnEl = nspx * nspy; % number of strain points per element

    % X direction first
    spInfo.spnEl = spnEl;
    spInfo.nstrain = nstrain;
    spInfo.nspx = nspx; spInfo.nspy = nspy;
end