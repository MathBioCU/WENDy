function w_cell = vec2cell(w,length_vec)
    ind=1;
    nstates = length(length_vec);
    w_cell = cell(nstates,1);
    for nn=1:nstates
        w_cell{nn} = w(ind:ind+length_vec(nn)-1);
        ind = ind+length_vec(nn);
    end
end
