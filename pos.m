function ix = pos(ix, Nx)
    if ix <= 0
       %ix=1;
       ix = Nx+ix;
    elseif ix > Nx
       ix=ix-Nx;
       %ix=Nx;
    end
end