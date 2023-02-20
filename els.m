
function w=els(G,b,batch_size,num_runs,meth) 
    [K,J] = size(G);
    W = zeros(J,num_runs);
    for rr=1:num_runs
        inds = randperm(K,floor(K/batch_size));
        W(:,rr) = G(inds,:) \ b(inds);
    end
    if isequal(meth,'median')
        w = median(W,2);
    elseif isequal(meth,'mean')
        w = mean(W,2);
    end
end