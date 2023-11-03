function u=TTR_denoise_func(x,th,Nx,Ny,Nz,tv_iters)

x=reshape(x,Nx,Ny,Nz);
u=tvdenoise(x,2/th,tv_iters);
u = u(:);