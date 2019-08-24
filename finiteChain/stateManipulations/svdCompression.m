function mps_out = svdCompression( mps_in, Dmax )
    %% Truncate bond dimnension of mps_in to Dmax using svd compression

    mps_out = rightLeftSweep_RCF(mps_in, Dmax);
end