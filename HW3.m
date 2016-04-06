%% Part 1 B&W
clear all; close all; clc

Abw=imread('derek2','jpeg');
[ny,nx]=size(Abw);
x0=double(uint8(nx/2)); y0=double(uint8(ny/2));
Abw=double(Abw);
kx=1:nx; ky=1:ny; 
[Kx,Ky]=meshgrid(kx,ky)
At=fft2(Abw);Ats=fftshift(At);

%Gaussian filter
fs=[0.01 0.001 0.0001 0];
figure(1)
for j=1:4
   F=exp(-fs(j)*(Kx-x0).^2-fs(j)*(Ky-y0).^2);
   Atsf=Ats.*F; Atf=ifftshift(Atsf); Af=ifft2(Atf);
   subplot(2,2,j), imshow(uint8(Af)), colormap(gray);
end

%% Part 2 RGB

A=imread('derek1','jpeg');
[ny,nx,nz]=size(A);
x0=double(uint8(nx/2)); y0=double(uint8(ny/2));
A=double(A);
kx=1:nx; ky=1:ny; 
[Kx,Ky]=meshgrid(kx,ky)

%Gaussian filter
fs=[0.01 0.001 0.0001 0];
figure(2);
for j=1:4
   F=exp(-fs(j)*(Kx-x0).^2-fs(j)*(Ky-y0).^2);
   for k=1:nz
      At=fft2(A(:,:,k));Ats=fftshift(At);
      Atsf(:,:,k)=Ats.*F; Atf(:,:,k)=ifftshift(Atsf(:,:,k)); Af(:,:,k)=ifft2(Atf(:,:,k));
   end
   subplot(2,2,j), imshow(uint8(Af));
end

%% Part 2 B&W

Bbw=imread('derek4','jpeg');
B2=Bbw(130:170,150:190);
B2=double(B2);
[nx,ny]=size(B2);
B2_2=reshape(B2,nx*ny,1);

x=linspace(0,1,nx); y=linspace(0,1,ny); 
dx=x(2)-x(1); dy=y(2)-y(1);
onex=ones(nx,1); oney=ones(ny,1);
Dx=(spdiags([onex -2*onex onex],[-1 0 1],nx,nx))/dx^2; Ix=eye(nx);
Dy=(spdiags([oney -2*oney oney],[-1 0 1],ny,ny))/dy^2; Iy=eye(ny);
L=kron(Iy,Dx)+kron(Dy,Ix);

tspan=[0 0.002 0.004 0.006]; D=0.2;
[t,Bsol]=ode113('image_rhs',tspan,B2_2,[],L,D);
for j=1:length(t)
   B2_clean=uint8(reshape(Bsol(j,:),nx,ny));
   Bbw(140:160,160:180)=B2_clean(11:31,11:31);
   subplot(2,2,j), imshow(Bbw);
end

%% Part 2 RGB

B=imread('derek3','jpeg');
B2=B(130:170,150:190,:);
B2=double(B2);
[nx,ny,nz]=size(B2);
for j=1:nz
   B2_2(:,:,j)=reshape(B2(:,:,j),nx*ny,1);
   
   x=linspace(0,1,nx); y=linspace(0,1,ny); 
   dx=x(2)-x(1); dy=y(2)-y(1);
   onex=ones(nx,1); oney=ones(ny,1);
   Dx=(spdiags([onex -2*onex onex],[-1 0 1],nx,nx))/dx^2; Ix=eye(nx);
   Dy=(spdiags([oney -2*oney oney],[-1 0 1],ny,ny))/dy^2; Iy=eye(ny);
   L=kron(Iy,Dx)+kron(Dy,Ix);
   
   tspan=[0 0.002 0.004 0.006]; D=0.2;
   [t,Bsol]=ode113('image_rhs',tspan,B2_2(:,:,j),[],L,D);
   sol(:,:,j)=Bsol;
end
for j=1:length(t)
   for k=1:nz
       Bsol=sol(:,:,k);
       B2_clean(:,:,k)=uint8(reshape(Bsol(j,:),nx,ny));
   end 
   B(140:160,160:180,:)=B2_clean(11:31,11:31,:);
   subplot(2,2,j), imshow(B);
end
