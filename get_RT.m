function [RT,L,Cov,lambda] = get_RT(L0,L1,w,diag_reg)
    if ~all(~w)
        L = L0+sum(L1.*reshape(w,1,1,[]),3);
    else
        L = L0;
    end
    Cov = L*L';
    RT = chol((1-diag_reg)*Cov+diag_reg*diag(diag(Cov)));
    lambda=diag_reg;
    RT = RT';
end



%     tau = 1/diag_reg;

%     [~,S,~]=svd(Cov,0);
%     S = diag(S);

%     S = sort(eigs(Cov,size(Cov,1)),'descend');
%     lambda = max((S(1)-tau*S(end))/(tau-1),0);
%     RT = chol(Cov+lambda*eye(size(Cov,1)));

%     lambda = min(diag(Cov))*diag_reg;
%     RT = chol(Cov+lambda*eye(size(Cov,1)));

%     lambda = max((S(1)-tau*S(end))/(tau*(min(diag(Cov))-S(end))-(max(diag(Cov))-S(1))),0);
%     RT = chol((1-lambda)*Cov+lambda*diag(diag(Cov)));
%     RT = RT/norm(RT);