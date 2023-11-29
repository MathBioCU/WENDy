param_length_vec = cellfun(@(x)length(x),features);
Jac_mat = build_Jac_sym(features,xobs);
[L0,L1] = get_Lfac(Jac_mat,param_length_vec,V_cell,Vp_cell);
dims = size(L1);
L = L0 + reshape(reshape(permute(L1,[3 1 2]),dims(3),[]).'*w_hat,dims(1),[]);