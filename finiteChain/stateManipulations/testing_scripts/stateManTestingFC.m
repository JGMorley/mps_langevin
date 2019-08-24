%% Main function to generate tests
% run by executing 'run(unitTesting);' in the Command Window
function tests = unitTesting
    tests = functiontests(localfunctions);
end

%% Test functions

function test_canonicalFormFC_LCF_1site(testCase)
    % random single-site input gets normalized correctly
    d = randi([2 12]);
    A_init = rand(d,1);
    norm_init = A_init'*A_init;
    
    mps = canonicalFormFC({0, {A_init, 1}}, 'LCF', true);
    A = mps{2}{1};
    
    norm = A'*A;
    verifyEqual(testCase,norm,1.,'AbsTol',1e-12)
    verifyEqual(testCase,A*sqrt(norm_init),A_init,'AbsTol',1e-12)
end
    
function test_canonicalFormFC_LCF_2sites(testCase)
    % random two-site input canonicalized correctly
    physDimList = randi([2 12], [1 2]);
    bondDimList = randi([2 12], [1 1]);
    mpsInit = {0, {rand(physDimList(1),bondDimList(1)), 1},...
              {rand(bondDimList(1),physDimList(2)),2}};
          
    mps = canonicalFormFC(mpsInit, 'LCF', true);
    L0 = mps{1}; [A1,L1] = mps{2}{:}; [A2,L2] = mps{3}{:};
 
    TOL = 1e-12;
    
    % state unchanged
    psi_before = mpsInit{2}{1}*mpsInit{3}{1};
    norm = trace(psi_before*psi_before');
    psi_before = psi_before / sqrt(norm);
    
    psi_after = A1*A2;
    verifyEqual(testCase,psi_before,psi_after,'AbsTol',TOL)
    
    % first site
    verifyEqual(testCase, L0,1.,'AbsTol',TOL)
    verifyEqual(testCase, trace(A1*L1*A1'),1.,'AbsTol',TOL)
    verifyEqual(testCase, A1'*A1,eye(size(A1,2)),'AbsTol',TOL)
    
    % second site
    verifyEqual(testCase, trace(A2*A2'),1., 'AbsTol', TOL)
    verifyEqual(testCase, A2*A2',L1, 'AbsTol',TOL)
    verifyEqual(testCase, L2,1.,'AbsTol',TOL)   
end

 function test_canonicalFormFC_LCF_8sites(testCase)
    % random 8-site input canonicalized correctly (right canonical form)
    physDimList = randi([2 12], [1 8]);
    bondDimList = randi([2 12], [1 7]);
    mpsInit = {0,...
           {rand(physDimList(1),bondDimList(1)), 1},...               % 1st site
           {rand(bondDimList(1),bondDimList(2),physDimList(2)),2},... % 2nd site
           {rand(bondDimList(2),bondDimList(3),physDimList(3)),3},... % 3rd site
           {rand(bondDimList(3),bondDimList(4),physDimList(4)),4},... % 4th site
           {rand(bondDimList(4),bondDimList(5),physDimList(5)),5},... % 5th site
           {rand(bondDimList(5),bondDimList(6),physDimList(6)),6},... % 6th site
           {rand(bondDimList(6),bondDimList(7),physDimList(7)),7},... % 7th site
           {rand(bondDimList(7),physDimList(8)),8}};                  % 8th site
    
    mps = canonicalFormFC(mpsInit, 'LCF', true);
    L0 = mps{1};         [A1,L1] = mps{2}{:}; [A2,L2] = mps{3}{:}; 
    [A3,L3] = mps{4}{:}; [A4,L4] = mps{5}{:}; [A5,L5] = mps{6}{:};    
    [A6,L6] = mps{7}{:}; [A7,L7] = mps{8}{:}; [A8,L8] = mps{9}{:};
    
    TOL = 1e-12;
    
    % Check density matrix unchanged
    A1I = mpsInit{2}{1}; A2I = mpsInit{3}{1}; A3I = mpsInit{4}{1};
    A4I = mpsInit{5}{1}; A5I = mpsInit{6}{1}; A6I = mpsInit{7}{1};
    A7I = mpsInit{8}{1}; A8I = mpsInit{9}{1};
    
    psi_before = ncon({A1I,A2I,A3I,A4I,A5I,A6I,A7I,A8I},...
                      {[-1 1],[1 2 -2],[2 3 -3],[3 4 -4],...
                       [4 5 -5],[5 6 -6],[6 7 -7],[7 -8]});
    norm = ncon({psi_before,conj(psi_before)}, {1:8,1:8});
    psi_before = psi_before / sqrt(norm);
    psi_after = ncon({A1,A2,A3,A4,A5,A6,A7,A8},...
                      {[-1 1],[1 2 -2],[2 3 -3],[3 4 -4],...
                       [4 5 -5],[5 6 -6],[6 7 -7],[7 -8]});
                   
    verifyEqual(testCase,psi_before,psi_after,'AbsTol',TOL)
                   
    % Check L0, L8
    verifyEqual(testCase, L0, 1., 'AbsTol', TOL)
    verifyEqual(testCase, L8, 1., 'AbsTol', TOL)
    
    % Check right conditions
    verifyEqual(testCase, A8*A8', L7, 'AbsTol', TOL)
    verifyEqual(testCase, ncon({A7, L7, conj(A7)},{[-1 1 2],[1 3],[-2 3 2]}),...
                          L6, 'AbsTol', TOL)
    verifyEqual(testCase, ncon({A6, L6, conj(A6)},{[-1 1 2],[1 3],[-2 3 2]}),...
                          L5, 'AbsTol', TOL)
    verifyEqual(testCase, ncon({A5, L5, conj(A5)},{[-1 1 2],[1 3],[-2 3 2]}),...
                          L4, 'AbsTol', TOL)
    verifyEqual(testCase, ncon({A4, L4, conj(A4)},{[-1 1 2],[1 3],[-2 3 2]}),...
                          L3, 'AbsTol', TOL)
    verifyEqual(testCase, ncon({A3, L3, conj(A3)},{[-1 1 2],[1 3],[-2 3 2]}),...
                          L2, 'AbsTol', TOL)
    verifyEqual(testCase, ncon({A2, L2, conj(A2)},{[-1 1 2],[1 3],[-2 3 2]}),...
                          L1, 'AbsTol', TOL)
    verifyEqual(testCase, trace(A1*L1*A1'), 1., 'AbsTol', TOL)
    
    % Check left conditions
    verifyEqual(testCase, A1'*A1, eye(size(A1,2)), 'AbsTol', TOL)
    verifyEqual(testCase, ncon({conj(A2),A2},{[1 -1 2],[1 -2 2]}),...
                          eye(size(A2,2)), 'AbsTol',TOL)
    verifyEqual(testCase, ncon({conj(A3),A3},{[1 -1 2],[1 -2 2]}),...
                          eye(size(A3,2)), 'AbsTol',TOL)
    verifyEqual(testCase, ncon({conj(A4),A4},{[1 -1 2],[1 -2 2]}),...
                          eye(size(A4,2)), 'AbsTol',TOL)
    verifyEqual(testCase, ncon({conj(A5),A5},{[1 -1 2],[1 -2 2]}),...
                          eye(size(A5,2)), 'AbsTol',TOL)
    verifyEqual(testCase, ncon({conj(A6),A6},{[1 -1 2],[1 -2 2]}),...
                          eye(size(A6,2)), 'AbsTol',TOL)
    verifyEqual(testCase, ncon({conj(A7),A7},{[1 -1 2],[1 -2 2]}),...
                          eye(size(A7,2)), 'AbsTol',TOL)
    verifyEqual(testCase, trace(A8*A8'), 1., 'AbsTol', TOL)
    
 end
 
function test_canonicalFormFC_RCF_1site(testCase)
    % random single-site input gets normalized correctly
    d = randi([2 12]);
    A_init = rand(d,1);
    norm_init = A_init'*A_init;
    
    mps = canonicalFormFC({0, {A_init, 1}});
    A = mps{2}{1};
    
    norm = A'*A;
    verifyEqual(testCase,norm,1.,'AbsTol',1e-12)
    verifyEqual(testCase,A*sqrt(norm_init),A_init,'AbsTol',1e-12)
end
    
function test_canonicalFormFC_RCF_2sites(testCase)
    % random two-site input canonicalized correctly
    physDimList = randi([2 12], [1 2]);
    bondDimList = randi([2 12], [1 1]);
    mpsInit = {0, {rand(physDimList(1),bondDimList(1)), 1},...
              {rand(bondDimList(1),physDimList(2)),2}};
          
    mps = canonicalFormFC(mpsInit);
    L0 = mps{1}; [A1,L1] = mps{2}{:}; [A2,L2] = mps{3}{:};
 
    TOL = 1e-12;
    
    % state unchanged
    psi_before = mpsInit{2}{1}*mpsInit{3}{1};
    norm = trace(psi_before*psi_before');
    psi_before = psi_before / sqrt(norm);
    
    psi_after = A1*A2;
    verifyEqual(testCase,psi_before,psi_after,'AbsTol',TOL)
    
    % first site
    verifyEqual(testCase, L0,1.,'AbsTol',TOL)
    verifyEqual(testCase, trace(A1*A1'),1.,'AbsTol',TOL)
    verifyEqual(testCase, A1'*A1,L1,'AbsTol',TOL)
    
    % second site
    verifyEqual(testCase, ncon({A2, conj(A2)}, {[-1 1 2], [-2 1 2]}),...
                          eye(size(A2,1)), 'AbsTol', TOL)
    verifyEqual(testCase, trace(A2'*L1*A2), 1., 'AbsTol',TOL)
    verifyEqual(testCase, L2,1.,'AbsTol',TOL)   
end

 function test_canonicalFormFC_RCF_8sites(testCase)
    % random 8-site input canonicalized correctly (right canonical form)
    physDimList = randi([2 12], [1 8]);
    bondDimList = randi([2 12], [1 7]);
    mpsInit = {0,...
           {rand(physDimList(1),bondDimList(1)), 1},...               % 1st site
           {rand(bondDimList(1),bondDimList(2),physDimList(2)),2},... % 2nd site
           {rand(bondDimList(2),bondDimList(3),physDimList(3)),3},... % 3rd site
           {rand(bondDimList(3),bondDimList(4),physDimList(4)),4},... % 4th site
           {rand(bondDimList(4),bondDimList(5),physDimList(5)),5},... % 5th site
           {rand(bondDimList(5),bondDimList(6),physDimList(6)),6},... % 6th site
           {rand(bondDimList(6),bondDimList(7),physDimList(7)),7},... % 7th site
           {rand(bondDimList(7),physDimList(8)),8}};                  % 8th site
    
    mps = canonicalFormFC(mpsInit);
    L0 = mps{1};         [A1,L1] = mps{2}{:}; [A2,L2] = mps{3}{:}; 
    [A3,L3] = mps{4}{:}; [A4,L4] = mps{5}{:}; [A5,L5] = mps{6}{:};    
    [A6,L6] = mps{7}{:}; [A7,L7] = mps{8}{:}; [A8,L8] = mps{9}{:};
    
    TOL = 1e-12;
    
    % Check density matrix unchanged
    A1I = mpsInit{2}{1}; A2I = mpsInit{3}{1}; A3I = mpsInit{4}{1};
    A4I = mpsInit{5}{1}; A5I = mpsInit{6}{1}; A6I = mpsInit{7}{1};
    A7I = mpsInit{8}{1}; A8I = mpsInit{9}{1};
    
    psi_before = ncon({A1I,A2I,A3I,A4I,A5I,A6I,A7I,A8I},...
                      {[-1 1],[1 2 -2],[2 3 -3],[3 4 -4],...
                       [4 5 -5],[5 6 -6],[6 7 -7],[7 -8]});
    norm = ncon({psi_before,conj(psi_before)}, {1:8,1:8});
    psi_before = psi_before / sqrt(norm);
    psi_after = ncon({A1,A2,A3,A4,A5,A6,A7,A8},...
                      {[-1 1],[1 2 -2],[2 3 -3],[3 4 -4],...
                       [4 5 -5],[5 6 -6],[6 7 -7],[7 -8]});
                   
    verifyEqual(testCase,psi_before,psi_after,'AbsTol',TOL)
                   
    % Check L0, L8
    verifyEqual(testCase, L0, 1., 'AbsTol', TOL)
    verifyEqual(testCase, L8, 1., 'AbsTol', TOL)
    
    % Check right conditions
    verifyEqual(testCase, A8*A8', eye(size(A8,1)), 'AbsTol', TOL)
    verifyEqual(testCase, ncon({A7, conj(A7)}, {[-1 1 2], [-2 1 2]}),...
                          eye(size(A7,1)), 'AbsTol', TOL)
    verifyEqual(testCase, ncon({A6, conj(A6)}, {[-1 1 2], [-2 1 2]}),...
                          eye(size(A6,1)), 'AbsTol', TOL)
    verifyEqual(testCase, ncon({A5, conj(A5)}, {[-1 1 2], [-2 1 2]}),...
                          eye(size(A5,1)), 'AbsTol', TOL)
    verifyEqual(testCase, ncon({A4, conj(A4)}, {[-1 1 2], [-2 1 2]}),...
                          eye(size(A4,1)), 'AbsTol', TOL)
    verifyEqual(testCase, ncon({A3, conj(A3)}, {[-1 1 2], [-2 1 2]}),...
                          eye(size(A3,1)), 'AbsTol', TOL)
    verifyEqual(testCase, ncon({A2, conj(A2)}, {[-1 1 2], [-2 1 2]}),...
                          eye(size(A2,1)), 'AbsTol', TOL)
    verifyEqual(testCase, trace(A1*A1'), 1., 'AbsTol', TOL)
    
    % Check left conditions
    verifyEqual(testCase, A1'*A1, L1, 'AbsTol', TOL)
    verifyEqual(testCase, ncon({conj(A2),L1,A2},{[1 -1 2], [3 1], [3 -2 2]}),...
                          L2, 'AbsTol',TOL)
    verifyEqual(testCase, ncon({conj(A3),L2,A3},{[1 -1 2], [3 1], [3 -2 2]}),...
                          L3, 'AbsTol',TOL)
    verifyEqual(testCase, ncon({conj(A4),L3,A4},{[1 -1 2], [3 1], [3 -2 2]}),...
                          L4, 'AbsTol',TOL)
    verifyEqual(testCase, ncon({conj(A5),L4,A5},{[1 -1 2], [3 1], [3 -2 2]}),...
                          L5, 'AbsTol',TOL)
    verifyEqual(testCase, ncon({conj(A6),L5,A6},{[1 -1 2], [3 1], [3 -2 2]}),...
                          L6, 'AbsTol',TOL)
    verifyEqual(testCase, ncon({conj(A7),L6,A7},{[1 -1 2], [3 1], [3 -2 2]}),...
                          L7, 'AbsTol',TOL)
    verifyEqual(testCase, trace(A8'*L7*A8), 1., 'AbsTol', TOL)
    
 end

function test_projectOutNullDims_unsortedError(testCase)
    % test for unsorted singular values error
    mps = {1, {[0 1;1 0], diag([0 1])}, {[1 0;0 1], 1}};
    
    verifyError(testCase, @()projectOutNullDims(mps),...
                'finiteChain:projectOutNullDims:notSorted');
end
    
function test_projectOutNullDims_nulldims(testCase)
    % fixture with null dims runs with no error and gives expected output
    B1 = .5*[1 1 0 0;1 -1 0 0];
    B2(:,:,1) = sqrt(2)*[0 1;1 0;0 0;0 0];
    B2(:,:,2) = sqrt(2)*[1 0;0 -1;0 0;0 0];
    B3(:,:,1) = .5*[1 1 0 0 0;-1 1 0 0 0];
    B3(:,:,2) = .5*[-1 1 0 0 0;-1 -1 0 0 0];
    B4 = [0 1;1 0;0 0;0 0;0 0];
    mps_withNullDims = {1,{B1,diag([.5,.5,0,0])},{B2,diag([.5,.5])},...
                          {B3,diag([.5,.5,0,0,0])},{B4,1}};
    
    mps = projectOutNullDims(mps_withNullDims);
    L0 = mps{1};         [A1,L1] = mps{2}{:}; [A2,L2] = mps{3}{:};
    [A3,L3] = mps{4}{:}; [A4,L4] = mps{5}{:};
    
    verifyEqual(testCase,L0,1)
    verifyEqual(testCase,A1,[.5 .5;.5 -.5]);
    verifyEqual(testCase,L1,diag([.5 .5]));
    verifyEqual(testCase,A2,cat(3,B2(1:2,:,1),B2(1:2,:,2)));
    verifyEqual(testCase,L2,diag([.5,.5]));
    verifyEqual(testCase,A3,cat(3,B3(:,1:2,1),B3(:,1:2,2)));
    verifyEqual(testCase,L3,diag([.5,.5]));
    verifyEqual(testCase,A4,[0 1;1 0]);
    verifyEqual(testCase,L4,1);
end

function test_rightLeftSweep_RCF(testCase)    
    % random input -> correct right condition satisfied
    % N=5 sites
    physDimList = randi([2 12], [1 5]);
    bondDimList = randi([2 12], [1 4]);
    mps = {0,...
           {rand(physDimList(1),bondDimList(1)), 1},...               % 1st site
           {rand(bondDimList(1),bondDimList(2),physDimList(2)),2},... % 2nd site
           {rand(bondDimList(2),bondDimList(3),physDimList(3)),3},... % 3rd site
           {rand(bondDimList(3),bondDimList(4),physDimList(4)),4},... % 4th site
           {rand(bondDimList(4),physDimList(5)),5}};                  % 5th site
   
    mps = rightLeftSweep_RCF(mps);
    A1 = mps{2}{1}; A2 = mps{3}{1}; A3 = mps{4}{1}; 
    A4 = mps{5}{1}; A5 = mps{6}{1};
    
    newBondDims = [size(A1,2),size(A2,2),size(A3,2),size(A4,2)];
    
    TOL = 1e-12;
    
    verifyEqual(testCase,trace(A1*A1'),1.,'AbsTol',TOL);
    verifyEqual(testCase,ncon({conj(A2),A2}, {[-1 1 2], [-2 1 2]}),...
                         eye(newBondDims(1)),'AbsTol',TOL);
    verifyEqual(testCase,ncon({conj(A3),A3}, {[-1 1 2], [-2 1 2]}),...
                         eye(newBondDims(2)),'AbsTol',TOL);
    verifyEqual(testCase,ncon({conj(A4),A4}, {[-1 1 2], [-2 1 2]}),...
                         eye(newBondDims(3)),'AbsTol',TOL);
    verifyEqual(testCase,A5*A5',eye(newBondDims(4)),'AbsTol',TOL);
end
    

function test_leftRightSweep_RCF_RCF(testCase)
    % random input -> correct left condition satisfied
    % N=5 sites
    physDimList = randi([2 12], [1 5]);
    bondDimList = randi([2 12], [1 4]);
    mps = {0,...
           {rand(physDimList(1),bondDimList(1)), 1},...               % 1st site
           {rand(bondDimList(1),bondDimList(2),physDimList(2)),2},... % 2nd site
           {rand(bondDimList(2),bondDimList(3),physDimList(3)),3},... % 3rd site
           {rand(bondDimList(3),bondDimList(4),physDimList(4)),4},... % 4th site
           {rand(bondDimList(4),physDimList(5)),5}};                  % 5th site
    
    mps = leftRightSweep_RCF(mps);
    L0 = mps{1};         [A1,L1] = mps{2}{:}; [A2,L2] = mps{3}{:}; 
    [A3,L3] = mps{4}{:}; [A4,L4] = mps{5}{:}; [A5,L5] = mps{6}{:};
    
    TOL = L5 * 10 * eps; % L5 is state norm
    
    verifyEqual(testCase,A1'*L0*A1,L1,'AbsTol',TOL);
    verifyEqual(testCase,ncon({conj(A2),L1,A2}, {[1 -1 3], [1 2], [2 -2 3]}),...
                         L2,'AbsTol',TOL);
    verifyEqual(testCase,ncon({conj(A3),L2,A3}, {[1 -1 3], [1 2], [2 -2 3]}),...
                         L3,'AbsTol',TOL);
    verifyEqual(testCase,ncon({conj(A4),L3,A4}, {[1 -1 3], [1 2], [2 -2 3]}),...
                         L4,'AbsTol',TOL);
    verifyEqual(testCase,trace(A5'*L4*A5),L5,'AbsTol',1e-10);
end

function test_checkValidMPS_errors(testCase)
    % check bad inputs raise the right errors
    
    % 1. Non cell-array input
    mps = {0, {rand(5,1), 0}};
    verifyError(testCase, @()checkValidMPS(cell2table(mps)), ...
                'finiteChain:checkValidMPS:notACell')
            
    % 2. Not a row cell-array
    mps = {0; {rand(5,1), 0}}; % NB the rogue semicolon!
    verifyError(testCase, @()checkValidMPS(mps), ...
                'finiteChain:checkValidMPS:not1row')
            
    % 3. Not of form {-,{-,-},...,{-,-}}
    mps = {0, {rand(3,2), 0}, {rand(2,2,3)}, {rand(2,3), 0}};
    verifyError(testCase, @()checkValidMPS(mps), ...
                'finiteChain:checkValidMPS:WrongShapeCellarray')
                
    % 4. Incompatible bond dimensions
    mps{3} = {rand(2,3,3), 0};
    verifyError(testCase, @()checkValidMPS(mps), ...
                'finiteChain:checkValidMPS:incompatibleBondDims')
end

function test_checkValidMPS(~)
    % check a few (inc 1-site) valid MPS pass ok
    
    physDimList = [2 2 3 2 2 2]; % so N=6 sites
    bondDimList = [2 3 2 4 2];
    mps = {0,...
           {rand(physDimList(1),bondDimList(1)), 1},...               % 1st site
           {rand(bondDimList(1),bondDimList(2),physDimList(2)),2},... % 2nd site
           {rand(bondDimList(2),bondDimList(3),physDimList(3)),3},... % 3rd site
           {rand(bondDimList(3),bondDimList(4),physDimList(4)),4},... % 4th site
           {rand(bondDimList(4),bondDimList(5),physDimList(5)),5},... % 5th site
           {rand(bondDimList(5),physDimList(6)),6}};                  % 6th site
    checkValidMPS(mps);
    
    mps = {0, {rand(5,1), 0}}; % 1 site, d = 5
    checkValidMPS(mps);
end

function test_mpsGeneratorFC_error(testCase)
    % check we get the right errors if we feed in bad inputs
    ds = [2 2 2];
    Ds = [2 2 2];
    
    verifyError(testCase, @()mpsGeneratorFC(ds,Ds), ...
                'mpsGeneratorFC:DsWrongLength')
end
    
function test_mpsGeneratorFC_singlesite(testCase)
    % check we get the right thing if we ask for a single-site MPS
    ds = 5;
    
    mps = mpsGeneratorFC(ds);
    
    verifyEqual(testCase, iscell(mps), logical(1)) % verify we have a cellarray
    verifyEqual(testCase, size(mps), [1 2])        % check size of 2nd element
    verifyEqual(testCase, size(mps{2}{1}), [5 1]); % check size of tensor
end

function test_mpsGeneratorFC_output(testCase)
    % check we get the right sort of object with the right sort of
    % properties  
    ds = [2 4 5 3];
    Ds = [12 8 2];
    
    mps = mpsGeneratorFC(ds, Ds);
    
    verifyEqual(testCase, iscell(mps), logical(1)) % verify we have a cellarray,
    verifyEqual(testCase, size(mps), [1 5]);       % ...of right size overall,
    verifyEqual(testCase, size(mps{1}), [1 1]);      % ...and for 1st element,
    for k=2:4
        verifyEqual(testCase, size(mps{k}), [1 2]);  % ...and for 2:N elements,
    end
    verifyEqual(testCase, size(mps{2}{1}), [2 12])   % ...and for 1st tensor,
    verifyEqual(testCase, size(mps{3}{1}), [12 8 4]) % ...and for 2nd tensor,
    verifyEqual(testCase, size(mps{4}{1}), [8 2 5])  % ...and for 3rd tensor,
    verifyEqual(testCase, size(mps{5}{1}), [2 3])    % ...and for 4th tensor.  
end        