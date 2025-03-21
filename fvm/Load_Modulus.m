function E_load = Load_Modulus(meshInfo,rf,F)
%%%===========================Copyright==================================%%%
	%%%   Version Nov. 2024
	%%%
	%%%   Gengxuan Zhu <zhugx@zju.edu.cn>
	%%%   Institute of Applied Mechanics,Zhejiang University
	%%%
	%%%===========================Description================================%%%
	%%% This is a function to calculate modulus on the load edge.
	%%%======================================================================%%%
    
    % Assemble
    Ivar = meshInfo.rightIvar;
    ndof = meshInfo.ndof;
    nNod = meshInfo.nNod;
    rf_C = fsparse(Ivar(:,1),Ivar(:,2),reshape(rf',[],1),[ndof,nNod]);
    nelx = meshInfo.nelx;
    nodx = nelx + 1;
    dof_iter = 2 * (nelx + 1);
    
    % load_rf = rf_C((dof_iter - 1):dof_iter:end - 1,nodx:nodx:end);
    % f = F((dof_iter - 1):dof_iter:end - 1);
    % E_load = load_rf \ f;
    

    nely = meshInfo.nely;
    nody = nely + 1;
    load_rf = sparse(2*nody,nody);
    f = sparse(2*nody,1);
    load_rf(1:2:end,:) = rf_C((dof_iter - 1):dof_iter:end,nodx:nodx:end);
    load_rf(2:2:end,:) = rf_C(dof_iter:dof_iter:end,nodx:nodx:end);
    f(1:2:end - 1,:) = F((dof_iter - 1):dof_iter:end - 1);
    f(2:2:end,:) = F(dof_iter:dof_iter:end);
    E_load = (load_rf' * load_rf) \ (load_rf' * f);
end