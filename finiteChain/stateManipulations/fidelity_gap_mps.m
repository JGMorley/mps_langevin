function fidelity_gap = fidelity_gap_mps( mps1, mps2 )
    %% calculate the fidelity gap between two mps
    %    || |1> - |2> || = <1|1> + <2|2> - 2*Re{<1|2>}
    
    fidelity_gap = fidelity_mps(mps1,mps1) + fidelity_mps(mps2,mps2);
    fidelity_gap = fidelity_gap - 2*real(fidelity_mps(mps1,mps2));
end