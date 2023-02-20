function phifun = dphi(phifun,diff_order)
    syms y; 
    phifun = matlabFunction(diff(phifun(y),diff_order));
end