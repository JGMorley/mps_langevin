function checkValidMPS( mpsIn, Dis1_ok )
    %% checkValidMPS Check input is a valid MPS cell-array, raise error if not

    if nargin==1, Dis1_ok=true; end
    
    %% 1. Check we have a cell array
    if ~iscell(mpsIn)
        err.identifier = 'finiteChain:checkValidMPS:notACell';
        err.message = 'Not a cell!';
        error(err)
    end

    %% 2. Check mpsIn is of the form {-,{-,-},...,{-,-}}
    size_mpsIn = size(mpsIn);
    oneRow = size_mpsIn(1) == 1;
    if ~oneRow
        err.identifier = 'finiteChain:checkValidMPS:not1row';
        err.message = 'MPS cell array must have only one row';
        error(err)
    end

    correctShape = 1;
    wrongShapeElements = [];
    for k = 2:size_mpsIn(2)
        if ~isequal(size(mpsIn{k}), [1 2])
            correctShape = 0;
            wrongShapeElements = [wrongShapeElements, num2str(k),', '];
        end
    end

    if ~correctShape
        err.identifier = 'finiteChain:checkValidMPS:WrongShapeCellarray';
        err.message = ['Elements ',wrongShapeElements(1:end-2),...
                       ' of MPS cell array have incorrect shape'];
        error(err);           
    end
    
    if isequal(size(mpsIn),[1 2]) 
        % N=1
        return
    end
    
    %% 3. Check neighbouring bond dimensions are compatible
    any_D_is_1 = 0;
    compatibleBondDims = 1;
    incompatiblePairs = [];
    dimk = size(mpsIn{1,2}{1},2); % dim of bond index of first tensor
    if dimk==1, any_D_is_1 = 1; end
    for k = 3:size_mpsIn(2) % loop over sites 2,...,N
        dimkplus = size(mpsIn{1,k}{1},1); % dim of first bond index
        if dimk ~= dimkplus
            compatibleBondDims = 0;
            pair = ['(',num2str(k-2),',',num2str(k-1),'), '];
            incompatiblePairs = [incompatiblePairs,pair];
        end
        dimk = size(mpsIn{1,k}{1},2);
        if dimk==1, any_D_is_1 = 1; end
    end
    
    if size_mpsIn(2) == 2 % treat single-site separately
        if dimk ~= 1
            compatibleBondDims = 0;
            pair = '(1,), ';
            incompatiblePairs = [incompatiblePairs,pair];
        end
    end

    if ~compatibleBondDims
        err.identifier = 'finiteChain:checkValidMPS:incompatibleBondDims';
        err.message = ['Pairs of sites ',incompatiblePairs(1:end-2),...
                       ' have incompatible bond dimensions'];
        error(err);
    end 
    
    if ~Dis1_ok && any_D_is_1
        err.identifier = 'finiteChain:checkValidMPS:Dis1_notok';
        err.message = 'At least one bond dimension = 1 and Dis1_ok=0';
        error(err);
    end
end

