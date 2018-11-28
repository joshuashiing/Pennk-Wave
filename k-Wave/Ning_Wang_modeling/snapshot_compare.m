clear;clc;

tmax=300;
nx=240;
nz=240;
dx=10;dz=10;

xx=0:dx:(nx-1)*dx;
zz=0:dz:(nz-1)*dz;


fid2=fopen('snap_chen.dat','r');
fid4=fopen('snap_xing.dat','r');
fid5=fopen('snap_zhu.dat','r');

  

for nt=1:tmax

    p2=fread(fid2,[nx nz],'float');
    p4=fread(fid4,[nx nz],'float');
    p5=fread(fid5,[nx nz],'float');

    pp1=(p2-p4);
    pp2=(p5-p4);
    pp3=(p5-p2);

end

    figure(55)
    imagesc(xx,zz,pp1)
    axis image
    h = gca; % ��ȡ��ǰ��ͼ�����ָ��
    set(h,'FontSize',16); % �������ִ�С��ͬʱӰ���������ע��ͼ��������ȡ�
    title('Difference between Chen & Xing');
    ylabel('Depth (m)'); % ������
    xlabel('Distance (m)'); % ������

    colorbar


    figure(66)
    imagesc(xx,zz,pp2)
    axis image
    h = gca; % ��ȡ��ǰ��ͼ�����ָ��
    set(h,'FontSize',16); % �������ִ�С��ͬʱӰ���������ע��ͼ��������ȡ�
    title('Difference between Zhu & Xing');
    ylabel('Depth (m)'); % ������
    xlabel('Distance (m)'); % ������

    colorbar
    
    figure(77)
    imagesc(xx,zz,pp3)
    axis image
    h = gca; % ��ȡ��ǰ��ͼ�����ָ��
    set(h,'FontSize',16); % �������ִ�С��ͬʱӰ���������ע��ͼ��������ȡ�
    title('Difference between Zhu & Chen');
    ylabel('Depth (m)'); % ������
    xlabel('Distance (m)'); % ������

    colorbar   
   
fclose(fid2);
fclose(fid4);
fclose(fid5);