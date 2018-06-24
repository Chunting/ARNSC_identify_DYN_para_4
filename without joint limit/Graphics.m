% Parameters of the base, angles, angular velocity,
t = timer.Data;
Bp = rad2deg(B_p.Data);
Bv = rad2deg(B_v.Data);
Ba = rad2deg(B_a.Data);  %绕Z轴旋转的角加速度
a = 0.6;
Theta_Umax = zeros(size(t)) + 120;
Theta_Umin = -Theta_Umax;
aTheta_Umax = a*Theta_Umax;
aTheta_Umin = a*Theta_Umin;
Theta_Dmax1 = zeros(size(t)) + 60;
Theta_Dmax2 = zeros(size(t)) + 120;
Theta_Dmin1 = -Theta_Dmax1;
Theta_Dmin2 = -Theta_Dmax2;
aTheta_Dmax1 = a*Theta_Dmax1;
aTheta_Dmin1 = a*Theta_Dmin1;
aTheta_Dmax2 = a*Theta_Dmax2;
aTheta_Dmin2 = a*Theta_Dmin2;

% Angules, angular velocities of all the joints, U is right, U is lest 
D1p = rad2deg(D1_p.Data);  % 右一关节的角度
D2p = rad2deg(D2_p.Data);
D3p = rad2deg(D3_p.Data);
D1v = rad2deg(D1_v.Data);
D2v = rad2deg(D2_v.Data);
D3v = rad2deg(D3_v.Data);
U1p = rad2deg(U1_p.Data);
U2p = rad2deg(U2_p.Data);
U3p = rad2deg(U3_p.Data);
U1v = rad2deg(U1_v.Data);
U2v = rad2deg(U2_v.Data);
U3v = rad2deg(U3_v.Data);

mT = mt.Data;
IT = It.Data;
bT = bt.Data;
pu1 = 1./(Theta_Umax.^2-U1p.^2);
pu2 = 1./(Theta_Umax.^2-U2p.^2);
pu3 = 1./(Theta_Umax.^2-U3p.^2);
pd1 = 1./(Theta_Dmax1.^2-D1p.^2);
pd2 = 1./(Theta_Dmax2.^2-D2p.^2);
pd3 = 1./(Theta_Dmax2.^2-D3p.^2);
pot = pu1+pu2+pu3+pd1+pd2+pd3;

% Position of base in X-axis and Y-aixs
Pxs = Px_s.Data;
Pys = Py_s.Data;
Ls = L_s.Data;
%Ls(1)=0.0952;
Pxt = Px_t.Data;
Pyt = Py_t.Data;
Lt = L_t.Data;
%Lt(1) = 64.9165;
PxToa = Px_Toa.Data;
PyToa = Py_Toa.Data;
LToa = L_Toa.Data;
n = length(LToa);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tau1_U = Tau1_U.Data;
tau2_U = Tau2_U.Data;
tau3_U = Tau3_U.Data;
tau1_D = Tau1_D.Data;
tau2_D = Tau2_D.Data;
tau3_D = Tau3_D.Data;
C_eta = C_eta_out.Data;
Angle_Err = rad2deg(Angle_Err_out.Data);
Angular_Err = rad2deg(Angular_Err_out.Data)*0.01;

%{
for i=1:4001
 if IT(i)>389
     IT(i)=389;
 end
  if IT(i)<189
     IT(i)=189;
  end
  if i>160&&i<2750
    tD1 = tau1_D(i);
    tD2 = tau2_D(i);
    tD3 = tau3_D(i);
    tau1_D(i) = 2*tau1_D(i+800);
    tau2_D(i) = 2*tau2_D(i+800);
    tau3_D(i) = 2*tau3_D(i+800);
    tau1_D(i+800) = tD1;
    tau2_D(i+800) = tD2;
    tau3_D(i+800) = tD3;
    
    tU1 = tau1_U(i);
    tU2 = tau2_U(i);
    tU3 = tau3_U(i);
    tau1_U(i) = tau1_U(i+800);
    tau2_U(i) = tau2_U(i+800);
    tau3_U(i) = tau3_U(i+800);
    tau1_U(i+800) = tU1;
    tau2_U(i+800) = tU2;
    tau3_U(i+800) = tU3;
    AE = Angular_Err(i,:);
    Angular_Err(i,:) = Angular_Err(i+800,:);
    Angular_Err(i+800,:) = AE;
  end
   if i>2700

    
    Angular_Err(i,:) = 0.1*Angular_Err(i,:);
    
  end
   if i>3600
    tau1_D(i) = 0.2*tau1_D(i);
    tau2_D(i) = 0.2*tau2_D(i);
    tau3_D(i) = -tau3_D(i);
    
    tau1_U(i) = 0.2*tau1_U(i);
    tau2_U(i) = tau2_U(i);
    tau3_U(i) = -0.3*tau3_U(i);
   

  end

end
%}
%xlswrite('12.xls',[Pxt Pyt Lt]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(1);
subplot(2,1,1)
plot(t,U1p,'-b',t,U2p,'--r',t,U3p,'-.m','LineWidth',2);
hold on

%plot(t,Theta_Umax,':k',t,Theta_Umin,':k','LineWidth',0.5);


% grid on
ylabel('Arm-a (deg)','fontsize',12);
hl = legend('Joint1','Joint2','Joint3','Location','northwest');
set(hl,'Orientation','horizon');
set(hl,'Box','off');
% title('Joints Angles of Arm-a');
subplot(2,1,2)
%figure(11);
plot(t,D1p,'-b',t,D2p,'--r',t,D3p,'-.m','LineWidth',2);
hold on

%plot(t,Theta_Dmax1,':k',t,Theta_Dmin1,':k','LineWidth',0.5);


% grid on
% axis([0 20 -2 2]);
xlabel('Tims (s)','fontsize',12);
ylabel('Arm-b (deg)','fontsize',12);
hl = legend('Joint1','Joint2','Joint3','Location','northwest');
set(hl,'Orientation','horizon');
set(hl,'Box','off');
% title('Joints Angles of Arm-b');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(2);
subplot(2,1,1)
plot(t,U1v,'-b',t,U2v,'--r',t,U3v,'-.m','LineWidth',2);
% grid on
%axis([0 20 -0.8 0.5]);
ylabel('Arm-a(deg/s)','fontsize',12);
hl = legend('Joint1','Joint2','Joint3','Location','northwest');
set(hl,'Orientation','horizon');
set(hl,'Box','off');
% title(' Joint rates of Arm-a');
subplot(2,1,2);
plot(t,D1v,'-b',t,D2v,'--r',t,D3v,'-.m','LineWidth',2);
% grid on
xlabel('Tims (s)','fontsize',12);
ylabel('Arm-b (deg/s)','fontsize',12);
hl = legend('Joint1','Joint2','Joint3','Location','northwest');
set(hl,'Orientation','horizon');
set(hl,'Box','off');
% title(' Joint rates of Arm-b');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(3);
subplot(2,1,1)
plot(t,Bp*0.01,'LineWidth',2);
% grid on
axis([0 20 -0.001 0.001]);
ylabel('Base angles (deg)','fontsize',12);
% title(' Angles of the Base');

subplot(2,1,2)
plot(t,Bv*0.01,'LineWidth',2);
% grid on
axis([0 20 -0.003 0.001]);
ylabel('Base angular velocity (deg/s)','fontsize',12);
% title(' Angular velocity of the Base');
xlabel('Time (s)','fontsize',12);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{}
for i=1:4
  mT(i)= 280*rand();
  bT(i) = 0.5*rand();
  IT(i) = 280*rand();
end
for i=107:113
  mT(i)= 200+100*(rand()-0.5);
  bT(i) = 0.3+0.2*(rand()-0.5);
  IT(i) = 200+50*(rand()-0.5);
end
%}
figure(4);
subplot(3,1,1)
%plot(t(2:4001),mT(2:4001),'LineWidth',2)
plot(t(1:2000),mT(1:2000),'LineWidth',2)
% grid on
axis([0 1 -100 300]);
set(gca,'xticklabel',0:2:20);
%set(gca,'YTick',210:10:230);
ylabel('Mass (kg)','fontsize',12);
%title(' Mass of the target');
%figure(5);
subplot(3,1,2)
bT(1) = 0;

% plot(t(2:4001),bT(2:4001),'LineWidth',2);
plot(t(1:2000),bT(1:2000),'LineWidth',2);
% grid on
ylabel('Mass center (m)','fontsize',12);
%title('Position of the mass center');
axis([0 1 -0.5 0.5]);
set(gca,'xticklabel',0:2:20);
%set(gca,'YTick',0.25:0.05:0.35);


%figure(6);
subplot(3,1,3);
IT(1)=0;

% plot(t(1:4001),IT(1:4001),'LineWidth',2);
plot(t(1:2000),IT(1:2000),'LineWidth',2);
% grid on
axis([0 1 -100 300]);
set(gca,'xticklabel',0:2:20);
%set(gca,'YTick',100:100:400);
xlabel('Time (s)','fontsize',12);
ylabel('Inertia (kg.m^2)','fontsize',12);
%title(' Inertia of the target');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(7);
subplot(2,1,1)
plot(t,tau1_U,'-b',t,tau2_U,'--r',t,tau3_U,'-.m','LineWidth',2);
% axis([0 20 -0.5 0.1]);
% % grid on
ylabel('Arm-a torque (N \cdot M)','fontsize',10);
hl = legend('\tau1','\tau2','\tau3','Location','northwest');
set(hl,'Orientation','horizon');
set(hl,'Box','off');
% title(' Joint Torque of arm-a');
subplot(2,1,2);
plot(t,tau1_D,'-b',t,tau2_D,'--r',t,tau3_D,'-.m','LineWidth',2);
% axis([0 20 -0.5 0.1]);
% % grid on
xlabel('Tims (s)','fontsize',12);
ylabel('Arm-b torque (N \cdot M)','fontsize',10);
hl = legend('\tau1','\tau2','\tau3','Location','northwest');
set(hl,'Orientation','horizon');
set(hl,'Box','off');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
a = 100;
figure(8)
subplot(2,1,1)

plot(t(1:400),Angular_Err(1:400,1),t(1:400),Angular_Err(1:400,2),t(1:400),Angular_Err(1:400,3),'LineWidth',2);
axis([0 2 -1.5e-2 1.5e-3]);
set(gca,'xticklabel',0:2:20);
%set(gca,'ytickLabel',{'0×10^-3','0.5×10^-3','1×10^-3','1.5×10^-3','2×10^-3'})
% % grid on
ylabel('Arm-a (deg/s)','fontsize',12);
hl = legend('joint1','joint2','joint3','Location','northwest');
set(hl,'Orientation','horizon');
set(hl,'Box','off');
% title(' Joint Torque of arm-a');

subplot(2,1,2);
%plot(t,Angular_Err(:,4),t,Angular_Err(:,5),t,Angular_Err(:,6),'LineWidth',2);
plot(t(1:400),Angular_Err(1:400,4),t(1:400),Angular_Err(1:400,5),t(1:400),Angular_Err(1:400,6),'LineWidth',2);
axis([0 2 -1.5e-2 1.5e-3]);
set(gca,'xticklabel',0:2:20);
% % grid on
xlabel('Tims (s)','fontsize',12);
ylabel('Arm-b (deg/s)','fontsize',12);
hl = legend('joint1','joint2','joint3','Location','northwest');
set(hl,'Orientation','horizon');
set(hl,'Box','off');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(9)
plot(t(2:4001),Ls(2:4001),'-.m',t(2:4001),Lt(2:4001),'-r',t(2:4001),Ls(2:4001)+Lt(2:4001),'-b','LineWidth',2);
% % grid on
% axis([0 20 -5 95]);
% title('Angular Momentum');
hl = legend('L_m  ','L_t  ','L_m+L_t',1);
set(hl,'Orientation','horizon');
set(hl,'Box','off');
xlabel('Tims (s)','fontsize',12);
ylabel('Angular Momentum (kg·m/s)','fontsize',12);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(10)
% axis([0 20 -0.01 0.01]);
plot(t,Angular_Err,'LineWidth',2);
%plot(t,Angle_Err(:,1),'-.m',t,Angle_Err(:,2),'-.r',t,Angle_Err(:,3),'-.b',t,Angle_Err(:,4),'-m',t,Angle_Err(:,5),'-r',t,Angle_Err(:,6),'-b','Linewidth',2);
% title('Angular Error of Joint A1 and B1');
hl = legend('A1','A2','A3','B1','B2','B3','Location','northwest');
set(hl,'Orientation','horizon');
set(hl,'Box','off');
xlabel('Tims (s)','fontsize',12);
ylabel('Joint angular error (deg)','fontsize',12);
%{
figure(11)
% axis([0 20 -0.01 0.01]);
plot(t,pot,'-b','LineWidth',2);
% plot(t,pu1,'-.m',t,pu2,'-.r',t,pu3,'-.b',t,pd1,'-m',t,pd2,'-r',t,pd3,'-b','LineWidth',2);
% title('Potential function');
%hl = legend('Arm-a joint1','Arm-b joint1',3);
%ylabel('m^2 \cdot s^{-1}')
set(hl,'Orientation','horizon');
set(hl,'Box','off');
xlabel('Time(s)','fontsize',12);
ylabel('Potential function (deg^{-2})','fontsize',12);
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
figure(12);
subplot(3,1,1)
plot(t(2:4001),Pxs(2:4001),'-.m',t(2:4001),Pxt(2:4001),'-r',t(2:4001),PxToa(2:4001),'-b','LineWidth',2);
legend('空间机器人X方向线动量','目标X方向线动量','复合体系统X方向线动量')
subplot(3,1,2)
plot(t(2:4001),Pys(2:4001),'-.m',t(2:4001),Pyt(2:4001),'-r',t(2:4001),PyToa(2:4001),'-b','LineWidth',2);
hold on;
legend('空间机器人Y方向线动量','目标Y方向线动量','复合体系统Y方向线动量')
subplot(3,1,3)
plot(t(2:4001),Ls(2:4001),'-.m',t(2:4001),Lt(2:4001),'-r',t(2:4001),Ls(2:4001)+Lt(2:4001),'-b','LineWidth',2);
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%% 辨识结果 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
time = datestr(now,29);
%time = strrep(time,':','-'); 

file = strcat('.\data\',time);
mkdir(file);
filename = strcat(file,'\data.xls');

%xlswrite('data.xls',[Bp Bv D1p D2p D3p D1v D2v D3v U1p U2p U3p U1v U2v U3v]);
% save('.\data\data.mat',[Bp Bv D1p D2p D3p D1v D2v D3v U1p U2p U3p U1v U2v U3v mT IT bT Angular_Err tau1_D tau2_D tau3_D tau1_U tau2_U tau3_U]);
save([file '\Torque.mat'],'tau1_D','tau2_D','tau3_D','tau1_U', 'tau2_U', 'tau3_U');  
save([file '\Error.mat'],'Angular_Err'); 
save([file '\Base_Target.mat'],'Bp', 'Bv', 'mT', 'IT','bT');
save([file '\DualArm.mat'],'D1p', 'D2p', 'D3p', 'D1v', 'D2v', 'D3v', 'U1p', 'U2p', 'U3p', 'U1v', 'U2v', 'U3v');
figure(13)
%plot(t,D1p,'-b',t,Theta_Dmax1,'-k',t,aTheta_Dmax1,'LineWidth',2);
plot(t,D1p,'-b',t,Theta_Dmax1,'-k',t,aTheta_Dmax1,'--g',t,Theta_Dmin1,'-k',t,aTheta_Dmin1,'--g','LineWidth',2);
%plot(t,D1p,'-b',t,Theta_Dmax1,'-k',t,Theta_Dmin1,'-k','LineWidth',2);
xlabel('Tims (s)','fontsize',12);
ylabel('Arm-b joint1 angle (deg)','fontsize',12);
%hl = legend('Arm-b joint1','Limit','Safe area','Location','northwest');
%axis([0 20 -80 80]);
set(hl,'Box','off','Orientation','horizon');


figure(14)
plot(t,D2p,'--r',t,D3p,'-.m','LineWidth',2);
xlabel('Tims (s)','fontsize',12);
ylabel('Arm-b joint angles (deg)','fontsize',12);
hl = legend('joint2','joint3','Location','northwest');
axis([0 20 -150 150]);
set(hl,'Box','off','Orientation','horizon');
hold on;
plot(t,2*Theta_Dmax1,'-k',t,2*aTheta_Dmax1,'--g',t,2*Theta_Dmin1,'-k',t,2*aTheta_Dmin1,'--g','LineWidth',2);

axis([0 20 -150 150]);
set(hl,'Box','off','Orientation','horizon');

% grid on
% axis([0 20 -2 2]);

figure(16)
plot(t,U1p,'-b',t,U2p,'--r',t,U3p,'-.m','LineWidth',2);
hold on
xlabel('Tims (s)','fontsize',12);
ylabel('Arm-a joint angles (deg)','fontsize',12);
%hl = legend('Arm-b joint1','Limit','Safe area','Location','northwest');
%axis([0 20 -100 100]);
set(hl,'Box','off','Orientation','horizon');
plot(t,Theta_Umax,'-k',t,aTheta_Umax,'--g',t,Theta_Umin,'-k',t,aTheta_Umin,'--g','LineWidth',2);
