function mpsOut = mpsGeneratorFC(  ds, Ds )
    %% mpsGeneratorFC Create random MPS by populating elements with rand()
    %    mpsOut = mpsGeneratorFC( spinDimList, bondDimList )
    %        Create MPS with spin dims and bond dims given by spinDimList
    %        and bondDimList respectively, and with complex elements populated 
    %        by rand()
    %%
    
    mpsOut{1} = 0;
    
    nSites = size(ds, 2);
    
    if nSites == 1
        mpsOut{2} = {rand([ds,1]), 0};
        return
    end
    
    if size(Ds, 2) ~= (nSites - 1)
        err.identifier = 'mpsGeneratorFC:DsWrongLength';
        err.message = 'array Ds has a number of elements not equal to nSites-1';
        error(err);
    end
    
    %% 1. First site
    dims = [ds(1), Ds(1)];
    mpsOut{2} = {rand(dims) + rand(dims)*1i, 0};
    
    if nSites == 2
       dims = [Ds(1), ds(2)];
       mpsOut{3} = {rand(dims) + rand(dims)*1i, 0};
       return
    end
    
    %% 2. Middle sites
    for n = 2:(nSites - 1)
        dims = [Ds(n-1), Ds(n), ds(n)];
        mpsOut{n + 1} = {rand(dims) + rand(dims)*1i, 0};
    end
    
    %% 3. Last site
    dims = [Ds(nSites - 1), ds(nSites)];
    mpsOut{nSites + 1} = {rand(dims) + rand(dims)*1i, 0};
end