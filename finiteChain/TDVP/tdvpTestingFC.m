%% Main function to generate tests
% run by executing 'run(unitTesting);' in the Command Window
function tests = tdvpTestingFC
    tests = functiontests(localfunctions);
end

%% Test functions

% function test_tangentSpace_cleanrun_edgesite(~) <<< write this!

function test_tangentSpace_cleanrun_notedgesite(~)
    %% Check tangentSpace runs without error given a valid input
    A = rand(randi([2 12], [1 3]));
    [D1,D2,~] = size(A);
    right = rand(D2) + 1i*rand(D2); right = (right + right')/2;
    left  = rand(D1) + 1i*rand(D1); left  = (left + left')/2;
    
    tangentSpace(2,A,left,right,'null');
    tangentSpace(2,A,left,right); % without method argument explicitly included
end

function test_tangentSpace_notfullrank(testCase)
    %% Check a not-full-rank input throws the right error
    A = rand(randi([2 12], [1 3])); % A[i,j,\sigma]
    [D1,D2,~] = size(A);
    right = eye(D2); left = eye(D1);
    
    % Duplicate one column of A[(i,\sigma), j] thereby reducing rank by one
    A(:,1,:) = A(:,2,:);
    
    verifyError(testCase, @()tangentSpace(2,A,left,right),...
                'finiteChain:tangentSpace:wrongNumberNullDims');
end

function test_calculateC_FC_fictures_Neq1(testCase)
    %% check we get the expected output for N=1 sites
    N = 1;
    d = randi([2 12]);
    A = rand([d 1]) + 1i*rand([d 1]);
    H = rand(d) + 1i*rand(d);
    mps = {0, {A,0}};
    
    verifyEqual(testCase,calculateC_FC(1,mps,H,N),(H*A).');
end

function test_calculateC_FC_fixtures_Ngt1(testCase)
    %% check legLinksList and tensorList against fixture results for N>1
    H = -1;
    
    % 2 sites
    N = 2; 
    mps = {0, {1i,0}, {2i,0}};
    
    n = 1; 
    [~,tensorList,legLinksList] = calculateC_FC(n,mps,H,N,'True');
    verifyEqual(testCase,tensorList,{1i,2i,-1,-2i});
    verifyEqual(testCase,legLinksList,{[1 2],[2 3],[-1 4 1 3],[-2 4]});
    
    n = 2;
    [~,tensorList,legLinksList] = calculateC_FC(n,mps,H,N,'True');
    verifyEqual(testCase,tensorList,{1i,2i,-1,-1i});
    verifyEqual(testCase,legLinksList,{[1 2],[2 3],[4 -2 1 3],[4 -1]});
    
    % 3 sites
    N = 3;
    mps = {0, {1i,0}, {2i,0}, {3i,0}};
    
    n = 1;
    [~,tensorList,legLinksList] = calculateC_FC(n,mps,H,N,'True');
    verifyEqual(testCase,tensorList,{1i,2i,3i,-1,-2i,-3i});
    verifyEqual(testCase,legLinksList,{[1 2],[2 4 3],[4 5],[-1 6 8 1 3 5],...
                                       [-2 7 6],[7 8]});
                                   
    n = 2;
    [~,tensorList,legLinksList] = calculateC_FC(n,mps,H,N,'True');
    verifyEqual(testCase,tensorList,{1i,2i,3i,-1,-1i,-3i});
    verifyEqual(testCase,legLinksList,{[1 2],[2 4 3],[4 5],[6 -3 7 1 3 5],...
                                       [6 -1],[-2 7]});
                                   
    n = 3;
    [~,tensorList,legLinksList] = calculateC_FC(n,mps,H,N,'True');
    verifyEqual(testCase,tensorList,{1i,2i,3i,-1,-1i,-2i});
    verifyEqual(testCase,legLinksList,{[1 2],[2 4 3],[4 5],[6 8 -2 1 3 5],...
                                       [6 7],[7 -1 8]});
    % 5 sites
    N = 5;
    mps = {0, {1i,0}, {2i,0}, {3i,0}, {4i,0}, {5i,0}};
    
    n = 1;
    [~,tensorList,legLinksList] = calculateC_FC(n,mps,H,N,'True');
    verifyEqual(testCase,tensorList,{1i,2i,3i,4i,5i,-1,-2i,-3i,-4i,-5i});
    verifyEqual(testCase,legLinksList,{[1 2],[2 4 3],[4 6 5],[6 8 7],[8 9],...
                                       [-1 10 12 14 16 1 3 5 7 9],[-2 11 10],...
                                       [11 13 12],[13 15 14],[15 16]});
                                   
    n = 3;
    [~,tensorList,legLinksList] = calculateC_FC(n,mps,H,N,'True');
    verifyEqual(testCase,tensorList,{1i,2i,3i,4i,5i,-1,-1i,-2i,-4i,-5i});
    verifyEqual(testCase,legLinksList,{[1 2],[2 4 3],[4 6 5],[6 8 7],[8 9],...
                                       [10 12 -3 13 15 1 3 5 7 9],[10 11],...
                                       [11 -1 12],[-2 14 13],[14 15]});
end