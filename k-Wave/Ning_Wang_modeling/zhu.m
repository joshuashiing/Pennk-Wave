%%%%%% 2D velocity-stress viscoacoustic wave modeling %%%%%
%%%%%% fractional laplacians operators %%%%%%%%
%%%% zhu's method %%%%%%%%%%
%%%%%%%%% homogeneous media %%%%%%%%%%

clear;clc;

dt=1e-3;%% time simple
fm=15; %%% domain frequency
tmax=350; %%time steps
to=1/fm;
t=-to:dt:(tmax-1)*dt-to;
rik=(1-2*(pi*fm.*t).^2).*exp(-(pi*fm.*t).^2);  %% wavelet

nxo=200;nzo=200;pml=20;
nx=nxo+2*pml;nz=nzo+2*pml;
dx=10;dz=10;
dkx=1/(nx*dx);dkz=1/(nz*dz);

vp=ones(nz,nx)*3000;  %velocity
dens=ones(nz,nx)*2000;%density
Qp=ones(nz,nx)*10; %% Q

xx=0:dx:(nx-1)*dx;
zz=0:dz:(nz-1)*dz;

kkx=zeros(nx,1);kkz=zeros(nz,1);
b1=1:1:nx;
b2=1:1:nz;

kkx(1:nx/2+1)=(2*pi*dkx*(b1(1:nx/2+1)-1));                
kkx(nx/2+2:end)=(2*pi*dkx*(nx+1-b1(nx/2+2:end)));
kkz(1:nz/2+1)=(2*pi*dkz*(b2(1:nz/2+1)-1));                 
kkz(nz/2+2:end)=(2*pi*dkz*(nz+1-b2(nz/2+2:end)));

k2w=zeros(nz,nx); % k^2
for j=1:nz
    for i=1:nx
        k2w(j,i)=(kkx(i)^2+kkz(j)^2);
    end
end
k2ww=k2w;
k2ww(1,1)=1e-8;

shiftx=ones(nx,1);
shiftz=ones(nz,1);
shiftx(1:nx/2+1)=exp(1i*pi/nx.*(b1(1:nx/2+1)-1));
shiftx(nx/2+2:nx)=exp(-1i*pi/nx.*(nx+1-b1(nx/2+2:nx)));
shiftz(1:nz/2+1)=exp(1i*pi/nz.*(b2(1:nz/2+1)-1));
shiftz(nz/2+2:nz)=exp(-1i*pi/nz.*(nz+1-b2(nz/2+2:nz)));

for i=pml+1:nxo+pml
   vp(1:pml,i)=vp(pml+1,i);
   vp(nzo+pml+1:nz,i)=vp(nzo+pml,i);
   
end 

for j=pml+1:nzo+pml
   vp(j,1:pml)=vp(j,pml+1);
   vp(j,nxo+pml+1:nx)=vp(j,nxo+pml);
end

vp(1:pml,1:pml)=vp(pml+1,pml+1);
vp(1:pml,nxo+pml+1:nx)=vp(pml+1,nxo+pml);
vp(nzo+pml+1:nz,1:pml)=vp(nzo+pml,pml+1);
vp(nzo+pml+1:nz,nxo+pml+1:nx)=vp(nzo+pml,nxo+pml);
     

bbx=ones(nz,nx);
bbz=ones(nz,nx);

for j=1:nz
    for i=1:nx-1
        bbx(j,i)=2/(dens(j,i+1)+dens(j,i));
    end
end
i=nx;
bbx(1:nz,i)=1./dens(1:nz,i);
for i=1:nx
    for j=1:nz-1
        bbz(j,i)=2/(dens(j+1,i)+dens(j,i));
    end
end
j=nz;
bbz(j,1:nx)=1./dens(j,1:nx);

rp=1/pi*atan(1./Qp);
wo=20*pi*fm;


Clam=dens.*vp.^2.*(cos(0.5*atan(1./Qp))).^2;
Clam=Clam./wo.^(2*rp);


yinda=Clam.*vp.^(2*rp).*cos(pi*rp);
tao=Clam.*vp.^(2*rp-1).*sin(pi*rp);



pw1=k2w.^rp;
pw2=k2ww.^(rp-0.5);


clear Clam Cmu Qp  vp  dens

p=zeros(nz,nx);
vx=zeros(nz,nx);vz=zeros(nz,nx);

midx1=zeros(nz,nx);midz1=zeros(nz,nx);
midx2=zeros(nz,nx);midz2=zeros(nz,nx);

midd2=zeros(nz,nx);
midxx1=zeros(nz,nx);
midzz1=zeros(nz,nx);
middd2=zeros(nz,nx);

zs=nz/2;xs=nx/2;

%%%%%%%%%% pml parameters %%%%%%%%%%
% kmax=0;
% alpha_max=pi*fm;                        %%% not equal zero value
% n1=2;n2=1;n3=1;
% R=1e-5;vmax=5000;
% widthx=pml*dx; 
% widthz=pml*dz;
% dmax_x=(1+n1+n2)*vmax*log(1/R)/(2*widthx);
% dmax_z=(1+n1+n2)*vmax*log(1/R)/(2*widthz);
% 
% int_dx_lf=zeros(1,pml);int_dx_rg=zeros(1,pml);
% half_dx_lf=zeros(1,pml);half_dx_rg=zeros(1,pml);
% 
% int_dz_up=zeros(1,pml);int_dz_dw=zeros(1,pml);
% half_dz_up=zeros(1,pml);half_dz_dw=zeros(1,pml);
% 
% int_kx_lf=zeros(1,pml);int_kx_rg=zeros(1,pml);
% half_kx_lf=zeros(1,pml);half_kx_rg=zeros(1,pml);
% 
% int_kz_up=zeros(1,pml);int_kz_dw=zeros(1,pml);
% half_kz_up=zeros(1,pml);half_kz_dw=zeros(1,pml);
% 
% int_alphax_lf=zeros(1,pml);int_alphax_rg=zeros(1,pml);
% half_alphax_lf=zeros(1,pml);half_alphax_rg=zeros(1,pml);
% 
% int_alphaz_up=zeros(1,pml);int_alphaz_dw=zeros(1,pml);
% half_alphaz_up=zeros(1,pml);half_alphaz_dw=zeros(1,pml);
% 
% %%%% int
% int_dx_lf(1:pml)=dmax_x.*(((pml-1)*dx-xx(1:pml))/widthx).^(n1+n2);
% int_dx_rg(1:pml)=dmax_x.*((xx(nx-pml+1:nx)-(nx-pml)*dx)/widthx).^(n1+n2);
% int_dz_up(1:pml)=dmax_z.*(((pml-1)*dz-zz(1:pml))/widthz).^(n1+n2);
% int_dz_dw(1:pml)=dmax_z.*((zz(nz-pml+1:nz)-(nz-pml)*dz)/widthz).^(n1+n2);
% 
% int_alphax_lf(1:pml)=alpha_max.*(xx(1:pml)./((pml-1)*dx)).^n3;
% int_alphax_rg(1:pml)=alpha_max.*(((nx-1)*dx-xx(nx-pml+1:nx))./((pml-1)*dx)).^n3;
% int_alphaz_up(1:pml)=alpha_max.*(zz(1:pml)./((pml-1)*dz)).^n3;
% int_alphaz_dw(1:pml)=alpha_max.*(((nz-1)*dz-zz(nz-pml+1:nz))./((pml-1)*dz)).^n3;
% 
% int_kx_lf(1:pml)=1+kmax.*(((pml-1)*dx-xx(1:pml))/widthx).^n1;
% int_kx_rg(1:pml)=1+kmax.*((xx(nx-pml+1:nx)-(nx-pml)*dx)/widthx).^n1;
% int_kz_up(1:pml)=1+kmax.*(((pml-1)*dz-zz(1:pml))/widthz).^n1;
% int_kz_dw(1:pml)=1+kmax.*((zz(nz-pml+1:nz)-(nz-pml)*dz)/widthz).^n1;
% %%%% half
% half_dx_lf(1:pml)=dmax_x.*((pml-0.5:-1:0.5)./pml).^(n1+n2);
% half_dx_rg(1:pml)=half_dx_lf(pml:-1:1);
% half_dz_up(1:pml)=dmax_z.*((pml-0.5:-1:0.5)./pml).^(n1+n2);
% half_dz_dw(1:pml)=half_dz_up(pml:-1:1);
% 
% half_alphax_lf(1:pml)=alpha_max.*((-0.5:1:pml-1.5)./(pml-1)).^n3;
% half_alphax_rg(1:pml)=half_alphax_lf(pml:-1:1);
% 
% half_alphaz_up(1:pml)=alpha_max.*((-0.5:1:pml-1.5)./(pml-1)).^n3;
% half_alphaz_dw(1:pml)=half_alphaz_up(pml:-1:1);
% 
% half_kx_lf(1:pml)=1+kmax.*((pml-0.5:-1:0.5)./pml).^n1;
% half_kx_rg(1:pml)=half_kx_lf(pml:-1:1);
% half_kz_up(1:pml)=1+kmax.*((pml-0.5:-1:0.5)./pml).^n1;
% half_kz_dw(1:pml)=half_kz_up(pml:-1:1);
% 
% fi_vx_lf=zeros(nz,pml);fi_vx_rg=zeros(nz,pml);
% fi_vz_up=zeros(pml,nx);fi_vz_dw=zeros(pml,nx);
% 
% fi_vx_up=zeros(pml,nx);fi_vx_dw=zeros(pml,nx);
% fi_vz_lf=zeros(nz,pml);fi_vz_rg=zeros(nz,pml);
% 
% fi_strxx_lf=zeros(nz,pml);fi_strxx_rg=zeros(nz,pml);
% fi_strzz_up=zeros(pml,nx);fi_strzz_dw=zeros(pml,nx);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

midxxzz_f=zeros(nz,nx);
midxz_f=zeros(nz,nx);
midx_f=zeros(nz,nx);
midz_f=zeros(nz,nx);
midzx_f=zeros(nz,nx);


recx=fopen('trace_zhu.dat','wb');

fx=fopen('snap_zhu.dat','wb');


for t=1:tmax
    
    %  vx  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fft_px=fft(p,[],2);    
   %%% rightward shift dx/2
    for j=1:nz
        midx1(j,1:nx)=((kkx./shiftx))'.*fft_px(j,1:nx);
        midx1(j,1:nx/2+1)=1i*midx1(j,1:nx/2+1);
        midx1(j,nx/2+2:nx)=-1i*midx1(j,nx/2+2:nx);
    end
    p_x=real(ifft(midx1,[],2));  % pxx to x

    
%      for j=1:nz
%         tem1=half_dx_lf(2:pml)./half_kx_lf(2:pml)+half_alphax_lf(2:pml); 
%         bxx=exp(-dt.*tem1);
%         tem2=(half_kx_lf(2:pml).^2).*tem1;
%         axx=half_dx_lf(2:pml)./tem2.*(bxx-1);
%         sumxx=p_x(j,1:pml-1);              %%%%
%         fi_strxx_lf(j,1:pml-1)=bxx.*fi_strxx_lf(j,1:pml-1)+axx.*sumxx;
%         tem3=(1-half_kx_lf(2:pml))./half_kx_lf(2:pml).*sumxx;
%         p_x(j,1:pml-1)=p_x(j,1:pml-1)+fi_strxx_lf(j,1:pml-1)+tem3;
%         
%         tem1=tem1(end:-1:1);
%         bxx=bxx(end:-1:1);
%         tem2=tem2(end:-1:1);
%         axx=axx(end:-1:1);
%         sumxx=p_x(j,nx-pml+1:nx-1);   %%%%
%         fi_strxx_rg(j,1:pml-1)=bxx.*fi_strxx_rg(j,1:pml-1)+axx.*sumxx;
%         tem3=(1-half_kx_rg(1:pml-1))./half_kx_rg(1:pml-1).*sumxx;
%         p_x(j,nx-pml+1:nx-1)=p_x(j,nx-pml+1:nx-1)+fi_strxx_rg(j,1:pml-1)+tem3;
%     end

    vx=vx+dt*bbx.*(p_x);

    
    %%%%%% vz %%%%%%%%%%%%%%%% 
    fft_pz=fft(p,[],1);   
    
    %%% upward shift dz/2
    for i=1:nx    
        midz2(1:nz,i)=((kkz.*shiftz)).*fft_pz(1:nz,i);  
        midz2(1:nz/2+1,i)=1i*midz2(1:nz/2+1,i);
        midz2(nz/2+2:nz,i)=-1i*midz2(nz/2+2:nz,i);
    end
    p_z=real(ifft(midz2,[],1)); % pzz to z
    

%     for i=1:nx
%         tem1=half_dz_up(2:pml)./half_kz_up(2:pml)+half_alphaz_up(2:pml);  
%         bzz=exp(-dt.*tem1);
%         tem2=(half_kz_up(2:pml).^2).*tem1;
%         azz=half_dz_up(2:pml)./tem2.*(bzz-1);
%         sumzz=p_z(1:pml-1,i);                                              
%         fi_strzz_up(1:pml-1,i)=(bzz').*fi_strzz_up(1:pml-1,i)+(azz').*sumzz;
%         tem3=(1-half_kz_up(2:pml))./half_kz_up(2:pml).*sumzz';
%         p_z(1:pml-1,i)=p_z(1:pml-1,i)+fi_strzz_up(1:pml-1,i)+tem3';
% 
%         
%         tem1=tem1(end:-1:1);
%         bzz=bzz(end:-1:1);
%         tem2=tem2(end:-1:1);
%         azz=azz(end:-1:1);
%         sumzz=p_z(nz-pml+1:nz-1,i);                              
%         fi_strzz_dw(1:pml-1,i)=(bzz').*fi_strzz_dw(1:pml-1,i)+(azz').*sumzz;
%         tem3=(1-half_kz_dw(1:pml-1))./half_kz_dw(1:pml-1).*sumzz';
%         p_z(nz-pml+1:nz-1,i)=p_z(nz-pml+1:nz-1,i)+fi_strzz_dw(1:pml-1,i)+tem3';
%     end


    vz=vz+dt*bbz.*(p_z);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% p   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
     fft_vxx=fft(vx,[],2);
     fft_vzz=fft(vz,[],1);
    
     %%% ���� dx/2
    for j=1:nz
        midx1(j,1:nx)=((kkx.*shiftx))'.*fft_vxx(j,1:nx);
        midx1(j,1:nx/2+1)=1i*midx1(j,1:nx/2+1);
        midx1(j,nx/2+2:nx)=-1i*midx1(j,nx/2+2:nx);
    end
     %%% ���� dz/2
    for i=1:nx    
        midz1(1:nz,i)=((kkz./shiftz)).*fft_vzz(1:nz,i);  
        midz1(1:nz/2+1,i)=1i*midz1(1:nz/2+1,i);
        midz1(nz/2+2:nz,i)=-1i*midz1(nz/2+2:nz,i);
    end
   
    mid1=real(ifft(midx1,[],2));  %%% vx to x
    mid2=real(ifft(midz1,[],1));  %%% vz to z
    
    %%%% calculate fi_vx and fi_vz
%      for i=1:nx
%         temp11=int_dz_up(1:pml)./int_kz_up(1:pml)+int_alphaz_up(1:pml);
%         bz=exp(-dt.*temp11);
%         temp22=(int_kz_up.^2).*temp11;
%         az=int_dz_up./temp22.*(bz-1);
%         sumz=mid2(1:pml,i);
%         fi_vz_up(1:pml,i)=(bz').*fi_vz_up(1:pml,i)+(az').*sumz;
%         temp33=(1-int_kz_up)./int_kz_up.*sumz';
%         mid2(1:pml,i)=mid2(1:pml,i)+fi_vz_up(1:pml,i)+temp33';
%         
%         temp11=temp11(end:-1:1);
%         bz=bz(end:-1:1);
%         temp22=temp22(end:-1:1);
%         az=az(end:-1:1);
%         sumz=mid2(nz-pml+1:nz,i);
%         fi_vz_dw(1:pml,i)=(bz').*fi_vz_dw(1:pml,i)+(az').*sumz;
%         temp33=(1-int_kz_dw)./int_kz_dw.*sumz';
%         mid2(nz-pml+1:nz,i)=mid2(nz-pml+1:nz,i)+fi_vz_dw(1:pml,i)+temp33';
%     end
%     
%     %%% disturb mid1, fi_vx_lf and fi_vx_rg
%      for j=1:nz
%         temp1=int_dx_lf(1:pml)./int_kx_lf(1:pml)+int_alphax_lf(1:pml);
%         bx=exp(-dt.*temp1);
%         temp2=(int_kx_lf.^2).*temp1;
%         ax=int_dx_lf./temp2.*(bx-1);
%         sumx=mid1(j,1:pml);
%         fi_vx_lf(j,1:pml)=bx.*fi_vx_lf(j,1:pml)+ax.*sumx;
%         temp3=(1-int_kx_lf)./int_kx_lf.*sumx;
%         mid1(j,1:pml)=mid1(j,1:pml)+fi_vx_lf(j,1:pml)+temp3;
%         
%         temp1=temp1(end:-1:1);
%         bx=bx(end:-1:1);
%         temp2=temp2(end:-1:1);
%         ax=ax(end:-1:1);
%         sumx=mid1(j,nx-pml+1:nx);
%         fi_vx_rg(j,1:pml)=bx.*fi_vx_rg(j,1:pml)+ax.*sumx;
%         temp3=(1-int_kx_rg)./int_kx_rg.*sumx;
%         mid1(j,nx-pml+1:nx)=mid1(j,nx-pml+1:nx)+fi_vx_rg(j,1:pml)+temp3;
%      end
    
     midxz=mid1+mid2;
     fft_mid3=fft2(midxz);
 
     md1p=yinda.*(real(ifft2(pw1.*fft_mid3)));
     md2p=tao.*(real(ifft2(pw2.*fft_mid3))-real(ifft2(pw2.*midxxzz_f)))/dt;
     
     midxxzz_f=fft_mid3;

     p=p+dt*(md1p+md2p);
     
     fwrite(fx,p,'float');
     
    p(zs,xs)=p(zs,xs)+rik(t);
    %%%%%%%%%%%%%%%%%%
    
    
    
 %%%%%%%%%   
%     mx=max(max(vx));
%     mn=min(min(vx));
    if(mod(t+5,5)==0)
        figure(18)

        imagesc(xx,zz,p)
        colorbar
        colormap(flipud(gray))
        axis image
%       caxis([mn mx]/1000)
        set(gcf,'color','w')
        title(['t=',num2str(t)],'fontsize',12)
        pause(0.001)
    end
       poux(t,:)=p(pml+80,pml+1:nx-pml);
      
    
 end
 fwrite(recx,poux,'float');


fclose(recx)
fclose(fx)
