function out=outn(f,in,n)
    s=cell(n,1);
    [s{:}] = f(in);
    out = s{n};
end