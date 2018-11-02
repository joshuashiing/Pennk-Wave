function [dx,dy,vp,Qh,rho]=padgridpml(velocity,Qtrue,rhob,r,npml,deltx,freesur)

% update history
% Add freesurface to make surface to be free to reflect
% Feb 22 2012 Tieyuan

%=========˥�������
    nxe=size(velocity,1);nze=size(velocity,2);
    x=1:nxe;
    z=1:nze;
    nn=npml;
    nx=nxe+2*nn;
    nz=nze+2*nn;
    %ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    %ccc     nn is the width grid of PML in the boundary
    %ccc     PML parameters
    %cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    deltan=nn*deltx; % Thickness of PML region
    d0=log(1/r)*3/2*max(max(velocity))/deltan;
    %---------PML attenuation function------------
    dx=zeros(nx,nz);dy=zeros(nx,nz);
    vp = zeros(nx,nz);Qh = zeros(nx,nz);rho = zeros(nx,nz);
    if strcmp(freesur,'no')
    % Extend velocity to include the PML region
    vp(nn+1:nxe+nn,nn+1:nze+nn)=velocity;
    Qh(nn+1:nxe+nn,nn+1:nze+nn)=Qtrue;
    rho(nn+1:nxe+nn,nn+1:nze+nn)=rhob;
    else
    % Extend velocity to include the PML region
    vp(1:nxe,nn+1:nze+nn)=velocity;
    Qh(1:nxe,nn+1:nze+nn)=Qtrue;
    rho(1:nxe,nn+1:nze+nn)=rhob;        
    end
    for  j=1:nn
        dy(:,j)=d0*((nn-j+1)/nn)^2;
        vp(:,j)=vp(:,nn+1);
        Qh(:,j)=Qh(:,nn+1);
        rho(:,j)=rho(:,nn+1);
    end

    for  j=(nz-nn):nz
        dy(:,j)=d0*((j-nz+nn)/nn)^2;
        vp(:,j)=vp(:,nz-nn-1);
        Qh(:,j)=Qh(:,nz-nn-1);
        rho(:,j)=rho(:,nz-nn-1);
    end
    
    if strcmp(freesur,'no')
    for i=1:nn
        dx(i,:)=d0*((nn-i+1)/nn)^2;
        %vp(i,:)=4.6*ones(size(vp(nn+1,:)));
        vp(i,:)=vp(nn+1+i-1,:);
        Qh(i,:)=Qh(nn+1+i-1,:);
        %rho(i,:)=2.2*ones(size(rho(nn+1,:)));
        rho(i,:)=rho(nn+1+i-1,:);
    end
    end
    k=0;
    for i=(nx-nn):nx
        k = k + 1;
        dx(i,:)=d0*((i-nx+nn)/nn)^2;
        %vp(i,:)=5.7*ones(size(vp(nx-nn-1,:)));
        vp(i,:)=vp(nx-nn-1,:);
        %Qh(i,:)=2.66*ones(size(Qh(nx-nn-1,:)));
        Qh(i,:)=Qh(nx-nn-1,:);
        %rho(i,:)=2.66*ones(size(rho(nx-nn-1,:)));
        rho(i,:)=rho(nx-nn-1,:);
    end
    

%         vp(nx-2*nn:nx,:)=vp(nx-4*nn:nx-2*nn,:);
%         Qh(nx-2*nn:nx,:)=Qh(nx-4*nn:nx-2*nn,:);
%         rho(nx-2*nn:nx,:)=rho(nx-4*nn:nx-2*nn,:);

 % number of PML boundary layers
% [vp,x3,z3] = padgrid(vp,x,z,2*npml+1);
% [Q,x3,z3] = padgrid(Q,x,z,2*npml+1);