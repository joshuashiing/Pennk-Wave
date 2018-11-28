clear;clc;

tmax=350;
dt=0.001;
tr=100;
nxo=200;nzo=200;pml=20;
nx=nxo+2*pml;nz=nzo+2*pml;
dx=10;dz=10;

xx=0:dx:(nx-1)*dx;
zz=0:dz:(nz-1)*dz;
tt=0:dt:(tmax-1)*dt;

fid1=fopen('trace_zhu.dat','r');
fid2=fopen('trace_xing.dat','r');
fid3=fopen('trace_chen.dat','r');


p1=fread(fid1,[tmax nx],'float');
p2=fread(fid2,[tmax nx],'float');
p3=fread(fid3,[tmax nx],'float');


trace1=p1(:,tr);
trace2=p2(:,tr);
trace3=p3(:,tr);

figure()
 
h = gca; % 
set(h,'FontSize',15); % 
hold on;
plot(tt,trace2,'k','LineWidth',1);
hold on;
plot(tt,trace1,'r','LineWidth',1);
hold on;
plot(tt,trace3,'b+','LineWidth',1);




ylabel('Amplitude'); % ï¿½ï¿½ï¿½ï¿½ï¿?
xlabel('Time (s)'); % ï¿½ï¿½ï¿½ï¿½ï¿?



 set(gcf,'Position',[100 350 860 269]);

   hg1 = legend('Xing', 'Zhu','Chen',0);

 
ah2 = axes('position',get(gca,'position'),'visible','off');  
