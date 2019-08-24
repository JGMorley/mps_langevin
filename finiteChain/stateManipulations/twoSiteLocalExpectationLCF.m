function [expn, lists] = twoSiteLocalExpectationLCF( O1, n1, O2, n2, mps_LCF )
    %% calculate expectation of two local operators on different sites:
    %  <O1,O2> where O1 acts on site n1 and O2 acts on site n2!=n1
    %
    %  lists has methods 'tensorList' and 'legLinksList', for debugging
  
    N = length(mps_LCF) - 1;
    
    if n1==n2
        error('n1=n2 not allowed')
    elseif n1>n2
        ntemp = n1; 
        n1 = n2;
        n2 = ntemp;
        
        Otemp = O1;
        O1 = O2;
        O2 = Otemp;
        
        clear ntemp Otemp
    end
    
    if n1==1 && n2<N
        nshift = n1 - 1; % subtract this to get indexing from 1 at n1-st site...
        n2Shft = n2 - nshift; % ...like so
        
        tensorList   = cell([1 2*n2Shft+3]);
        legLinksList = cell([1 2*n2Shft+3]);
        
        % build tensorList
        for n=n1:n2
            An = mps_LCF{n+1}{1};
            tensorList{n-nshift} = conj(An);
        end
        
        tensorList{n2Shft+1} = mps_LCF{n2+1}{2}; 
        tensorList{n2Shft+2} = O1; 
        tensorList{n2Shft+3} = O2;
        
        for n=n1:n2
            An = mps_LCF{n+1}{1};
            tensorList{n2Shft+n-nshift+3} = An;
        end
        
        % build legLinksList
        legLinksList{1}        = [1 2];                   % bottom-left
        legLinksList{n2Shft+1} = [3*n2Shft+1 2*n2Shft];   % Rn2
        legLinksList{n2Shft+2} = [1 2*n2Shft+1];          % O1
        legLinksList{n2Shft+3} = [2*n2Shft-1 3*n2Shft+2]; % O2
        legLinksList{n2Shft+4} = 2*n2Shft + [1 2];        % top-left
        
        for n=(n1+1):n2
            nShft = n - nshift;
            legLinksList{nShft} = 2*nShft - [2 0 1];
            legLinksList{n2Shft+3+nShft} = [2*n2Shft+nShft,  ...
                                            2*n2Shft+nShft+1,...
                                            2*nShft-1];
        end
        legLinksList{2*n2Shft+3} = 3*n2Shft + [0 1 2]; % overwrite top-right
        
        warning('off','ncon:suboptimalsequence')
        expn = ncon(tensorList,legLinksList);
        warning('on','ncon:suboptimalsequence')
        
        if nargout>1
            lists.tensorList = tensorList;
            lists.legLinksList = legLinksList;
        end
    elseif n1==1 && n2==N
        nshift = n1 - 1; % subtract this to get indexing from 1 at n1-st site...
        n2Shft = n2 - nshift; % ...like so
        
        tensorList   = cell([1 2*n2Shft+2]);
        legLinksList = cell([1 2*n2Shft+2]);
        
        % build tensorList
        for n=n1:n2
            An = mps_LCF{n+1}{1};
            tensorList{n-nshift} = conj(An);
        end
        
        tensorList{n2Shft+1} = O1; 
        tensorList{n2Shft+2} = O2;
        
        for n=n1:n2
            An = mps_LCF{n+1}{1};
            tensorList{n2Shft+n-nshift+2} = An;
        end
        
        % build legLinksList
        legLinksList{1}          = [1 2];                 % bottom-left
        legLinksList{n2Shft}     = 2*n2Shft - [2 1];      % bottom-right 
        legLinksList{n2Shft+1}   = [1 2*n2Shft];          % O1
        legLinksList{n2Shft+2}   = [2*n2Shft-1 3*n2Shft]; % O2
        legLinksList{n2Shft+3}   = 2*n2Shft + [0 1];      % top-left
        legLinksList{2*n2Shft+2} = 3*n2Shft - [1 0];      % top-right
        
        for n=(n1+1):(n2-1)
            nShft = n - nshift;
            legLinksList{nShft} = 2*nShft - [2 0 1];
            legLinksList{n2Shft+2+nShft} = [2*n2Shft+nShft-1,...
                                            2*n2Shft+nShft,  ...
                                            2*nShft-1];
        end
        
        warning('off','ncon:suboptimalsequence')
        expn = ncon(tensorList,legLinksList);
        warning('on','ncon:suboptimalsequence')
        
        if nargout>1
            lists.tensorList = tensorList;
            lists.legLinksList = legLinksList;
        end
    elseif n1>1 && n2<N
        nshift = n1 - 1; % subtract this to get indexing from 1 at n1-st site...
        n2Shft = n2 - nshift; % ...like so
        
        tensorList   = cell([1 2*n2Shft+3]);
        legLinksList = cell([1 2*n2Shft+3]);
        
        % build tensorList
        for n=n1:n2
            An = mps_LCF{n+1}{1};
            tensorList{n-nshift} = conj(An);
        end
        
        tensorList{n2Shft+1} = mps_LCF{n2+1}{2}; 
        tensorList{n2Shft+2} = O1; 
        tensorList{n2Shft+3} = O2;
        
        for n=n1:n2
            An = mps_LCF{n+1}{1};
            tensorList{n2Shft+n-nshift+3} = An;
        end
        
        % build legLinksList
        legLinksList{1}        = [2*n2Shft+1 2 1];          % bottom-left
        legLinksList{n2Shft+1} = [3*n2Shft+2 2*n2Shft];   % Rn2
        legLinksList{n2Shft+2} = [1 2*n2Shft+2];          % O1
        legLinksList{n2Shft+3} = [2*n2Shft-1 3*n2Shft+3]; % O2
        legLinksList{n2Shft+4} = 2*n2Shft + [1 3 2];      % top-left
        
        for n=(n1+1):n2
            nShft = n - nshift;
            legLinksList{nShft} = 2*nShft - [2 0 1];
            legLinksList{n2Shft+3+nShft} = [2*n2Shft+nShft+1,  ...
                                            2*n2Shft+nShft+2,...
                                            2*nShft-1];
        end
        legLinksList{2*n2Shft+3} = 3*n2Shft + [1 2 3]; % overwrite top-right
        
        warning('off','ncon:suboptimalsequence')
        expn = ncon(tensorList,legLinksList);
        warning('on','ncon:suboptimalsequence')
        
        if nargout>1
            lists.tensorList = tensorList;
            lists.legLinksList = legLinksList;
        end
    elseif n1>1 && n2==N
        nshift = n1 - 1; % subtract this to get indexing from 1 at n1-st site...
        n2Shft = n2 - nshift; % ...like so
        
        tensorList   = cell([1 2*n2Shft+2]);
        legLinksList = cell([1 2*n2Shft+2]);
        
        % build tensorList
        for n=n1:n2
            An = mps_LCF{n+1}{1};
            tensorList{n-nshift} = conj(An);
        end
        
        tensorList{n2Shft+1} = O1; 
        tensorList{n2Shft+2} = O2;
        
        for n=n1:n2
            An = mps_LCF{n+1}{1};
            tensorList{n2Shft+n-nshift+2} = An;
        end
        
        % build legLinksList
        legLinksList{1}          = [2*n2Shft 2 1];          % bottom-left
        legLinksList{n2Shft}     = 2*n2Shft - [2 1];        % bottom-right 
        legLinksList{n2Shft+1}   = [1 2*n2Shft+1];          % O1
        legLinksList{n2Shft+2}   = [2*n2Shft-1 3*n2Shft+1]; % O2
        legLinksList{n2Shft+3}   = 2*n2Shft + [0 2 1];      % top-left
        legLinksList{2*n2Shft+2} = 3*n2Shft + [0 1];        % top-right
        
        for n=(n1+1):(n2-1)
            nShft = n - nshift;
            legLinksList{nShft} = 2*nShft - [2 0 1];
            legLinksList{n2Shft+2+nShft} = [2*n2Shft+nShft,...
                                            2*n2Shft+nShft+1,  ...
                                            2*nShft-1];
        end
        
        warning('off','ncon:suboptimalsequence')
        expn = ncon(tensorList,legLinksList);
        warning('on','ncon:suboptimalsequence')
        
        if nargout>1
            lists.tensorList = tensorList;
            lists.legLinksList = legLinksList;
        end
    end
end

