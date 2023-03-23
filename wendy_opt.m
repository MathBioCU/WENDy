%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% WENDy: covariance-corrected ODE parameter estimation
%%%%%%%%%%%% Copyright 2023, All Rights Reserved
%%%%%%%%%%%% Code by Daniel Ames Messenger

%%% 'meth' options: 'LS', 'elsLS'
%%% 'LS' options: (none)
%%% 'elsLS' options: batch_size, num_runs, avg_meth

function w = windy_opt(G,b,varargin)

    defaultmeth = 'LS';
    defaultbatch_size = 1;
    defaultnum_runs = 1;
    defaultavg_meth = 'mean';
    defaultcov = [];
    
    inp = inputParser;
    addParameter(inp,'meth',defaultmeth);
    addParameter(inp,'batch_size',defaultbatch_size);
    addParameter(inp,'num_runs',defaultnum_runs);
    addParameter(inp,'avg_meth',defaultavg_meth);
    addParameter(inp,'cov',defaultcov);
    
    parse(inp,varargin{:});  
    
    meth = inp.Results.meth;
    batch_size = inp.Results.batch_size;
    num_runs = inp.Results.num_runs;
    avg_meth = inp.Results.avg_meth;
    cov = inp.Results.cov;

    if isequal(meth,'LS')
        if isempty(cov)
            w = G \ b;
        else
            w = lscov(G,b,cov);
        end
    end

    if isequal(meth,'TLS')
        [~,~,V] = svd([G b],0);
        [~,n]=size(G);
        w = V(1:n,n+1:end);
    end

    if isequal(meth,'ensLS')
        w = els(G,b,batch_size,num_runs,avg_meth);
    end

end
