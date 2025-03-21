function spValue = Sp_Value(meshInfo,spInfo,N,dNds,dNdt,locJ,eps)
%%%===========================Copyright==================================%%%
	%%%   Version Nov. 2024
	%%%
	%%%   Gengxuan Zhu <zhugx@zju.edu.cn>
	%%%   Institute of Applied Mechanics,Zhejiang University
	%%%
	%%%===========================Description================================%%%
	%%% This is a function to calculate the value of integral expression.
	%%%======================================================================%%%
    % Generate 8x4 cell of shape functions
    % to calculate 32 numerical integration     
    nEl = meshInfo.nEl;
    spnEl = spInfo.spnEl;
    N1 = repmat(N{1},8,1);
    N2 = repmat(N{2},8,1);
    N3 = repmat(N{3},8,1);
    N4 = repmat(N{4},8,1);
    cellN = mat2cell([N1,N2,N3,N4],nEl * ones(8,1),spnEl * ones(1,4));
    % Generate 8x4 cell of strain values
    % to calculate 32 numerical integration
    % Calculate c1,c3,c5,c7(along X direction)
    C1 = cellfun(@(dNds,dNdt) ...
         eps{1} .* (locJ{1} .* dNds + locJ{2} .* dNdt) +  ...
         eps{2} .* (locJ{3} .* dNds + locJ{4} .* dNdt), ...
         dNds,dNdt,'UniformOutput',false);
    c1 = repmat(C1{1},1,4);
    c3 = repmat(C1{2},1,4);
    c5 = repmat(C1{3},1,4);
    c7 = repmat(C1{4},1,4);
    % Calculate c2,c4,c6,c8(along Y direction)
    C2 = cellfun(@(dNds,dNdt) ...
         eps{3} .* (locJ{1} .* dNds + locJ{2} .* dNdt) +  ...
         eps{4} .* (locJ{3} .* dNds + locJ{4} .* dNdt), ...
         dNds,dNdt,'UniformOutput',false);
    c2 = repmat(C2{1},1,4);
    c4 = repmat(C2{2},1,4);
    c6 = repmat(C2{3},1,4);
    c8 = repmat(C2{4},1,4);
    cellC = mat2cell([c1;c2;c3;c4;c5;c6;c7;c8],nEl * ones(8,1),spnEl * ones(1,4));
    % Calculate integration values
    spValue = cellfun(@(cellN,cellC) cellN .* cellC,cellN,cellC,'UniformOutput',false);
end