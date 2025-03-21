function value = Double_Integral(spInfo,spValue)
%%%===========================Copyright==================================%%%
	%%%   Version Oct. 2024
	%%%
	%%%   Gengxuan Zhu <zhugx@zju.edu.cn>
	%%%   Institute of Applied Mechanics,Zhejiang University
	%%%
	%%%===========================Description================================%%%
	%%% This is a function to calculate double integral by trapezoidal rule.
	%%%======================================================================%%%
    nspx = spInfo.nspx;
    nspy = spInfo.nspy;
    spnEl = spInfo.spnEl;
    spsEl = spInfo.spsEl(:,1:nspx);
    sptEl = spInfo.sptEl(:,1:nspx:end);
    nEl = size(spsEl,1);
    h = (spsEl(:,2:end) - spsEl(:,1:end - 1)) / 2;
    k = (sptEl(:,2:end) - sptEl(:,1:end - 1)) / 2;
    % Coefficient Generation
    h_index(:,1) = h(:,1); h_index(:,nspx) = h(:,end);
    h_index(:,2:(end - 1)) = h(:,1:(end - 1)) + h(:,2:end);
    h_index = num2cell(h_index,2);
    k_index(:,1) = k(:,1); k_index(:,nspy) = k(:,end);
    k_index(:,2:(end - 1)) = k(:,1:(end - 1)) + k(:,2:end);
    k_index = num2cell(k_index,2);
    % Coefficient Matrix Generation
    H = cellfun(@(h_index) repmat(h_index',1,size(sptEl,2)),h_index,'UniformOutput',false);
    K = cellfun(@(k_index) repmat(k_index,size(spsEl,2),1),k_index,'UniformOutput',false);
    M = cell2mat(cellfun(@(H,K) reshape(H .* K,1,[]),H,K,'UniformOutput',false));
    M_tot = mat2cell(repmat(M,8,4),nEl * ones(1,8),spnEl * ones(1,4));
    % Calculate integration for each element
    value = cellfun(@(M_tot,spValue) sum(M_tot .* spValue,2),M_tot,spValue,'UniformOutput',false);
    % Convergence Analysis
    % xdata = linspace(-1,1,nspx);
    % ydata = linspace(-1,1,nspy);
    % [plot_x,plot_y] = meshgrid(xdata,ydata);
    % testdata = plot_x .^4 .* plot_y .^4;
    % test = mat2cell(repmat(repmat(reshape(testdata',1,[]),meshInfo.nEl,1),8,4), ...
    %                meshInfo.nEl * ones(1,8),spInfo.spnEl * ones(1,4));
end