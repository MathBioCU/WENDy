function dx = rhs_fun(features,params,x)
    nstates = length(x);
    x = num2cell(x);
    dx = zeros(nstates,1);
    for i=1:nstates
        dx(i) = cellfun(@(z1) z1(x{:}),features{i})*params{i}(:);
    end
end