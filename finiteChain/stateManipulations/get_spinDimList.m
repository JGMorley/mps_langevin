function spinDimList = get_spinDimList(mps)
    N = length(mps) - 1;
    spinDimList=zeros([1 N]);
    
    % n=1
    spinDimList(1) = size(mps{2}{1},1);
    
    % 1<n<N
    for n=2:(N-1)
        spinDimList(n) = size(mps{n+1}{1},3);
    end
    
    % n=N
    spinDimList(N) = size(mps{N+1}{1},2);
end

