function V = solNorm(U,S)
    MU = pagemtimes(S.MassMatrix,S.J.*U);
    
    MU1 = MU(:,:,1); U1 = U(:,:,1);
    MU2 = MU(:,:,2); U2 = U(:,:,2);
    MU3 = MU(:,:,3); U3 = U(:,:,3);
    
    L21 = diag(U1'*MU1); L22 = diag(U2'*MU2); L23 = diag (U3'*MU3);
    
    V = sqrt(sum(L21)+sum(L22)+sum(L23));
end