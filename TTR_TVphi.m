function y=TTR_TVphi(X,Nx,Ny,Nz)

X=reshape(X,Nx,Ny,Nz);

[y,dif]=TTR_TVnorm(X);

function [y,dif]=TTR_TVnorm(x)

TV=TTR_TV3D_conv(x);

dif=sqrt(sum(TV.*conj(TV),4));

y=sum(dif(:));
end
end