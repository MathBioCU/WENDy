%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% WENDy: covariance-corrected ODE parameter estimation
%%%%%%%%%%%% Copyright 2023, All Rights Reserved
%%%%%%%%%%%% Code by Daniel Ames Messenger

function [w_hat,res,res_0,w_hat_its,V_cell,Vp_cell,...
    Theta_cell,mt,xobs,Jac_mat,G_0,b_0,RT,stdW,mseW,CovW] = wendy_fcn_0(...
    xobs,tobs,features,toggle_smooth,...
    mt_cell,mt_min,mt_max,K_min,K_max,center_scheme,toggle_VVp_svd,...
    w0,optim_params,iter_diff_tol,max_iter,diag_reg,pvalmin,check_pval_it)

    %%% get dimensions
    param_length_vec = cellfun(@(x)length(x),features);
    eq_inds = cellfun(@(x) ~isempty(x), features);
    num_eq = sum(eq_inds);
    [M,nstates] = size(xobs);

    %%% estimate noise variance and build initial covariance 
    sig_ests = arrayfunvec(xobs,@estimate_sigma,1);
    RT_0 = spdiags(kron(sig_ests(:),ones(M,1)),0,M*nstates,M*nstates);

    %%% smooth data
    if toggle_smooth > 0
        init_m_fac = []; expand_fac = 1.5; max_filter_fac = []; maxits_SMAF = []; verbose = []; 
        sws = zeros(nstates,1);
        RT_0_temp=cell(nstates,1);
        for nn=1:nstates
            if toggle_smooth == 1
                [sws(nn),~,~] = get_optimal_SMAF(tobs,xobs(:,nn),[],init_m_fac,max_filter_fac,expand_fac,maxits_SMAF,[],verbose);
            else
                sws(nn) = toggle_smooth;
            end
            xobs(:,nn) = movmean(xobs(:,nn),2*sws(nn)+1,1,'Endpoints','shrink');
            RT_0_temp{nn} = spdiags(ones(M,2*sws(nn)+1),[-sws(nn):sws(nn)],M,M);
            RT_0_temp{nn} = RT_0_temp{nn}./sum(RT_0_temp{nn},2);
        end
        RT_0 = blkdiag(RT_0_temp{:})*RT_0;
    end

    %%% get test function
    [cm,cn] = size(mt_cell);
    K = min(floor(K_max/nstates/cm), M);

    if cn<nstates
        mt = cell2mat(cellfun(@(mtc)arrayfunvec(xobs,@(x)get_rad(x,tobs,mtc{1},mtc{2},mtc{3},mt_min,mt_max),1),mt_cell(:,1),'uni',0));
        mt = ceil(1./mean(1./mt,2));
        if and(toggle_VVp_svd~=0,cm>1)
            V = cell2mat(cellfun(@(x,y)get_VVp_svd(y,tobs,K,x{1},center_scheme),mt_cell,num2cell(mt),'uni',0));
            [V,Vp] = VVp_svd(V,K_min,tobs,toggle_VVp_svd);
        else
            [V,Vp] = cellfun(@(x,y)get_VVp(y,tobs,1,K,x{1},center_scheme),mt_cell,num2cell(mt),'uni',0);
            V = cat(1,V{:}); Vp = cat(1,Vp{:});
        end
        V_cell = repmat({V},nstates,1);
        Vp_cell = repmat({Vp},nstates,1);
        mt = mt*ones(1,nstates);
    else
        mt = cell2mat(cellfun(@(mtc,x)get_rad(x,tobs,mtc{1},mtc{2},mtc{3},mt_min,mt_max),mt_cell,repmat(mat2cell(xobs,M,ones(1,nstates)),cm,1),'uni',0));
        V_cell = cell(num_eq,1);
        Vp_cell = cell(num_eq,1);
        for nn=1:num_eq
            V_cell{nn} = [];
            Vp_cell{nn} = [];
            for j=1:length(find(mt(:,nn)))
                if and(toggle_VVp_svd~=0,cm>1)
                    V_cell{nn} = [V_cell{nn};get_VVp_svd(mt(j,nn),tobs,K,mt_cell{j,min(nn,end)}{1},center_scheme)];
                else
                    [V,Vp] = get_VVp(mt(j,nn),tobs,1,K,mt_cell{j,min(nn,end)}{1},center_scheme);
                    V_cell{nn} = [V_cell{nn};V];
                    Vp_cell{nn} = [Vp_cell{nn};Vp];
                end
            end
            if and(toggle_VVp_svd~=0,cm>1)
                [V_cell{nn},Vp_cell{nn}] = VVp_svd(V_cell{nn},K_min,tobs,toggle_VVp_svd);
            end
        end
    end

    %%% build linear system
    xobs_cell = mat2cell(xobs,M,ones(1,nstates));
    Theta_cell = cellfun(@(x) cell2mat(cellfun(@(y) y(xobs_cell{:}), x, 'uni',0)), features,'uni',0);
    G_0 = cellfun(@(x,y) x*y,V_cell(:),Theta_cell(:),'uni',0);
    G_0 = blkdiag(G_0{:});
    b_0 = cell2mat(cellfun(@(x,y) -x*y,Vp_cell(:),xobs_cell(:),'uni',0));
    
    %%% build library Jacobian
    Jac_mat = build_Jac_sym(features,xobs);
    if max_iter>1
        [L0,L1] = get_Lfac(Jac_mat,param_length_vec,V_cell,Vp_cell);
        L0=L0*RT_0;
        L1=pagemtimes(L1,full(RT_0));
    end    

    % initialize
    if isempty(w0)
        w0 = wendy_opt(G_0,b_0,optim_params{:});
    end
    w_hat = w0;
    w_hat_its = w_hat;
    res = G_0*w_hat-b_0;
    res_0 = res;
    iter = 1; check = 1; pval=1;
    RT = eye(size(b_0,1));
    [~,pvals,~]=swtest(res);

    while all([check>iter_diff_tol,iter<max_iter,pval>pvalmin])

        %%% update covariance
        [RT,~,~,~] = get_RT(L0,L1,w_hat,diag_reg);
        G = RT \ G_0;
        b = RT \ b_0;

        %%% update parameters
        w_hat = wendy_opt(G,b,optim_params{:});
        res_n = G*w_hat-b;

        %%% check stopping conditions
        [~,pvals(iter+1),~]=swtest(res_n);
        if iter+1>check_pval_it
            pval = pvals(iter+1);
        end
        check = norm(w_hat_its(:,end)-w_hat)/norm(w_hat_its(:,end));
        iter = iter + 1;
        
        %%% collect quantities of interest
        res = [res res_n];
        res_0 = [res_0 (G_0*w_hat-b_0)];
        w_hat_its = [w_hat_its w_hat];

    end

    if pval<pvalmin
        disp(['error: WENDy iterates diverged'])
        [~,ind] = max(pvals);
        w_hat = w_hat_its(:,ind);
        if ind==1
            RT = eye(size(b_0,1));
        elseif ind<size(w_hat_its,2)
            [RT,~,~,~] = get_RT(L0,L1,w_hat_its(:,ind-1),diag_reg);
        end
        res = [res res(:,ind)];
        res_0 = [res_0 res_0(:,ind)];
        w_hat_its = [w_hat_its w_hat_its(:,ind)];
    end

    Ginv = G_0 \ RT;
    CovW = Ginv*Ginv';
    stdW = sqrt(diag(CovW)); 
    mseW = rms(res(:,end))^2;

end