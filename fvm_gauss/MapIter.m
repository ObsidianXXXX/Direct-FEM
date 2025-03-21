function zeta = MapIter(vec,noddata,eps,iter_num)
%%%===========================Copyright==================================%%%
	%%%   Version Oct. 2024
	%%%
	%%%   Gengxuan Zhu <zhugx@zju.edu.cn>
	%%%   Institute of Applied Mechanics,Zhejiang University
	%%%
	%%%===========================Description================================%%%
	%%% This is a function to map strain points from global coordinates to
    %%% local coordinates by interations.
	%%%======================================================================%%%
    % Initial solution
    zeta0 = [0.5;-0.5];
    zeta_his = zeta0;
    % Calculate global point
    [N,dNds,dNdt] = Point_Shape_Fun(zeta0);
    vec0 = noddata * N';
    vec_his = vec0;
    for i = 1:iter_num
        J = [noddata * dNds',noddata * dNdt'];
        invJ = inv(J);
        zeta = zeta0 + invJ * (vec - vec0);
        zeta_his = [zeta_his,zeta];
        [N,dNds,dNdt] = Point_Shape_Fun(zeta);
        vec0 = noddata * N';
        vec_his = [vec_his,vec0];
        if norm(vec0 - vec) < eps
            break
        end
        zeta0 = zeta;
    end
end