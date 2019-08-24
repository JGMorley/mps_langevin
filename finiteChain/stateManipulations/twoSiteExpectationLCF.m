function [expn, lists] = twoSiteExpectationLCF( Onnp, n, mps_LCF )
    %% calculate expectation of a two-site operator
    %  <O_{n,n+1}> where O_{n,n+1} acts on site n and n+1
    %
    %  lists has methods 'tensorList' and 'legLinksList', for debugging

    % 1. extract variables and check dims of O
    N = length(mps_LCF) - 1;
    assert(n<N, sprintf('require 1 <= n < N: N=%i and n=%i',N,n))
    
    An  = mps_LCF{n+1}{1};
    Anp = mps_LCF{n+2}{1};
    Rnp = mps_LCF{n+2}{2};
    
    spinDimList = get_spinDimList(mps_LCF);
    dn  = spinDimList(n);
    dnp = spinDimList(n+1);
    
    if ~isequal(size(Onnp), [dn dnp dn dnp])
        errmsg = ['Operator on ',num2str(n),'th site has wrong dimensions:'];
        errmsg = [errmsg,newline,'size(O) = [',num2str(size(O)),'],',newline,...
                  'Local dimensions = ',num2str([dn dnp]),'.'];
        error(errmsg)
    end
    
    %% 2. create tensorList and legLinksList for ncon
    % four cases, depending on whether n, (n+1) are at end of chain
    
    if (n+1) == N
        tensorList = {An,Anp,conj(An),conj(Anp),Onnp};
        if n==1
            legLinksList = {[1 2],[2 3],[4 5],[5 6],[4 6 1 3]}; 
        else
            legLinksList = {[1 2 3],[2 4],[1 6 5],[6 7],[5 7 3 4]};
        end
    else
        tensorList = {An,Anp,Rnp,conj(An),conj(Anp),Onnp};
        if n==1
            legLinksList = {[1 3],[3 4 2],[4 5],[8 7],[7 5 6],[8 6 1 2]};
        else
            legLinksList = {[1 2 3],[2 8 5],[8 6],[1 7 4],[7 6 9],[4 9 3 5]};
        end
    end
    
    % 2. Calculate <O>
    expn = ncon(tensorList,legLinksList);
    
    if nargout==2
        lists.tensorList = tensorList;
        lists.legLinksList = legLinksList;
    end
    
    % 3. Raise warning if expectation has an imaginary part
    if abs(imag(expn)) > 10*eps
        warning('Expecatation has an imaginary part')
    end
end

