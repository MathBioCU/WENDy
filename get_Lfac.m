function [L0,L1] = get_Lfac(Jac_mat,Js,V_cell,Vp_cell)
    [~,d,M] = size(Jac_mat);
    Jac_mat = permute(Jac_mat,[2 3 1]);
    eq_inds = find(Js);
    num_eq = length(eq_inds);
    L0 = blkdiag(Vp_cell{:});
    L1 = repmat(L0*0,1,1,sum(Js));
    Ktot = 0;
    Jtot = 0;
    for i=1:num_eq
        [K,~] = size(V_cell{i});
        J = Js(eq_inds(i));
        for ell=1:d
            L1(Ktot+(1:K),(ell-1)*M+(1:M),Jtot+(1:J)) = Jac_mat(ell,:,Jtot+(1:J)).*V_cell{i};
        end
        Ktot = Ktot+K;
        Jtot = Jtot+J;
    end
end