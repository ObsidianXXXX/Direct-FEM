function [spsEl,sptEl] = Global2Local(meshInfo,spInfo,eps,iter_num)
%%%===========================Copyright==================================%%%
	%%%   Version Oct. 2024
	%%%
	%%%   Gengxuan Zhu <zhugx@zju.edu.cn>
	%%%   Institute of Applied Mechanics,Zhejiang University
	%%%
	%%%===========================Description================================%%%
	%%% This is a function to map strain points from global coordinates to
    %%% local coordinates.
	%%%======================================================================%%%
    spnEl = spInfo.spnEl;
    % Node Info
    xEle = meshInfo.xEle;
    yEle = meshInfo.yEle;
    xdata = repmat(xEle,1,spnEl);
    ydata = repmat(yEle,1,spnEl);
    % Strain point info
    spxEl = spInfo.spxEl;
    spyEl = spInfo.spyEl;
    % Element data
    vec = num2cell([reshape(spxEl',1,[]); ...
                    reshape(spyEl',1,[])],1);
    noddata = mat2cell([reshape(xdata',1,[]);reshape(ydata',1,[])], ...
                        2, ...
                        4 * ones(1,numel(spxEl)));
    % Get mapping result
    zeta = cellfun(@(vec,noddata) MapIter(vec,noddata,eps,iter_num), ...
                       vec,noddata,'UniformOutput',false);
    zeta_mat = cell2mat(zeta);
    sps = reshape(zeta_mat(1,:),spnEl,[]);
    spsEl = sps';
    spt = reshape(zeta_mat(2,:),spnEl,[]);
    sptEl = spt';
end
