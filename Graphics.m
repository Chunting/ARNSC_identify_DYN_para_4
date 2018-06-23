% Parameters of the base, angles, angular velocity,
t = timer.Data;
Bp = rad2deg(B_p.Data);
Bv = rad2deg(B_v.Data);
Ba = rad2deg(B_a.Data);  %绕Z轴旋转的角加速度

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
Angular_Err = rad2deg(Angular_Err_out.Data);
Lamda = Lamda_out.Data;


%%%%%%%%%%%%%%%%%%%%%   1. 关节角度  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1);
subplot(2,1,1)
plot(t,U1p,'-b',t,U2p,'--r',t,U3p,'-.m','LineWidth',2);
hold on
% grid on
%xlabel('Time (s)','fontsize',12);
ylabel('Arm-a (deg)','fontsize',12);
hl = legend('Joint1','Joint2','Joint3','Location','northwest');
set(hl,'Orientation','horizon');
set(hl,'Box','off');
% title('Joints Angles of Arm-a');
subplot(2,1,2)
%figure(2);
plot(t,D1p,'-b',t,D2p,'--r',t,D3p,'-.m','LineWidth',2);
hold on
% grid on
% axis([0 20 -2 2]);
xlabel('Time (s)','fontsize',12);
ylabel('Arm-b (deg)','fontsize',12);
hl = legend('Joint1','Joint2','Joint3','Location','northwest');
set(hl,'Orientation','horizon');
set(hl,'Box','off');
% title('Joints Angles of Arm-b');
%%%%%%%%%%%%%%%%%%%%%%%%%  2. 关节速度  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(3);
subplot(2,1,1)
plot(t,U1v,'-b',t,U2v,'--r',t,U3v,'-.m','LineWidth',2);
% grid on
%axis([0 20 -0.8 0.5]);
%xlabel('Time (s)','fontsize',12);
ylabel('Arm-a (deg/s)','fontsize',12);
hl = legend('Joint1','Joint2','Joint3','Location','northwest');
set(hl,'Orientation','horizon');
set(hl,'Box','off');
% title(' Joint rates of Arm-a');
subplot(2,1,2);
%figure(4);
plot(t,D1v,'-b',t,D2v,'--r',t,D3v,'-.m','LineWidth',2);
% grid on
xlabel('Time (s)','fontsize',12);
ylabel('Arm-b (deg/s)','fontsize',12);
hl = legend('Joint1','Joint2','Joint3','Location','northwest');
set(hl,'Orientation','horizon');
set(hl,'Box','off');
% title(' Joint rates of Arm-b');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  3. 基座扰动  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(5);
subplot(2,1,1)
plot(t,Bp,'LineWidth',2);
% grid on
%axis([0 20 -0.001 0.001]);
xlabel('Time (s)','fontsize',12);
ylabel('Base angles (deg)','fontsize',12);
% title(' Angles of the Base');
% Bv(1)=0;
subplot(2,1,2)
plot(t,Bv,'LineWidth',2);
% grid on
axis([0 20 -inf inf]);
ylabel('Base angular velocity (deg/s)','fontsize',12);
% title(' Angular velocity of the Base');
xlabel('Time (s)','fontsize',12);
%%%%%%%%%%%%%%%%%%%%%%%%%%   4.辨识结果   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(6);
subplot(3,1,1)

semilogx(t,mT,'LineWidth',2)
axis([0 20 100 150]);
ylabel('Mass (kg)','fontsize',12);

subplot(3,1,2)
semilogx(t,bT,'LineWidth',2);
ylabel('Mass center (m)','fontsize',12);
axis([0 20 -0.5 0.5]);

subplot(3,1,3);
semilogx(t,IT,'LineWidth',2);
axis([0 20 -inf inf]);
xlabel('Time (s)','fontsize',12);
ylabel('Inertia (kg.m^2)','fontsize',12);

%%%%%%%%%%%%%%%%%%%%%%   5. 关节力矩     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(7);
%subplot(2,1,1)
semilogx(t,tau1_U,'-b',t,tau2_U,'--r',t,tau3_U,'-.m','LineWidth',2);
axis([0 20 -inf inf]);
set(gca,'xticklabel',0:2:20);
xlabel('Time (s)','fontsize',12);
ylabel('Arm-a torque (N \cdot M)','fontsize',12);
hl = legend('\tau1','\tau2','\tau3','Location','northwest');
set(hl,'Orientation','horizon');
set(hl,'Box','off');
% title(' Joint Torque of arm-a');
%subplot(2,1,2);
figure(8);
semilogx(t,tau1_D,'-b',t,tau2_D,'--r',t,tau3_D,'-.m','LineWidth',2);
axis([0 20 -inf inf]);
set(gca,'xticklabel',0:2:20);
xlabel('Time (s)','fontsize',12);
ylabel('Arm-b torque (N \cdot M)','fontsize',12);
hl = legend('\tau1','\tau2','\tau3','Location','northwest');
set(hl,'Orientation','horizon');
set(hl,'Box','off');
%%%%%%%%%%%%%%%%%%%%%    6. 动量守恒    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(9)
l = length(t);
plot(t(2:l),Ls(2:l),'-.m',t(2:l),Lt(2:l),'-r',t(2:l),Ls(2:l)+Lt(2:l),'-b','LineWidth',2);
hl = legend('L_m  ','L_t  ','L_m+L_t');
set(hl,'Orientation','horizon');
set(hl,'Box','off');
xlabel('Time (s)','fontsize',12);
ylabel('Angular Momentum (kg·m/s)','fontsize',12);

%%%%%%%%%%%%%%%%%%%%%%  7. 关节角速度误差  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
a = 100;
figure(10);
%subplot(2,1,1)
semilogx(t,Angular_Err(:,1),t,Angular_Err(:,2),t,Angular_Err(:,3),'LineWidth',2);
axis([0 20 -inf inf]);
%set(gca,'xticklabel',0:2:20);
xlabel('Time (s)','fontsize',12);
ylabel('Arm-a','fontsize',12);
hl = legend('joint1','joint2','joint3','Location','northwest');
set(hl,'Orientation','horizon');
set(hl,'Box','off');

figure(11);
%semilogx(t,bT,'LineWidth',2);
semilogx(t,Angular_Err(:,4:6),'LineWidth',2);
axis([0 20 -inf inf]);
axis('auto y');
%set(gca,'xticklabel',0:2:20);
xlabel('Time (s)','fontsize',12);
ylabel('Arm-b ','fontsize',12);
hl = legend('joint1','joint2','joint3','Location','northwest');
set(hl,'Orientation','horizon');
set(hl,'Box','off');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(12)
semilogx(t,Angle_Err(:,1),'-.m',t,Angle_Err(:,2),'-.r',t,Angle_Err(:,3),'-.b',t,Angle_Err(:,4),'-m',t,Angle_Err(:,5),'-r',t,Angle_Err(:,6),'-b','Linewidth',2);
axis([0 20 -inf inf]);
% title('Angular Error of Joint A1 and B1');
hl = legend('A1','A2','A3','B1','B2','B3','Location','northwest');
set(hl,'Orientation','horizon');
set(hl,'Box','off');
xlabel('Time (s)','fontsize',12);
ylabel('Joint angular error (deg)','fontsize',12);

%%%%%%%%%%%%%%%%%%%  9. 线动量  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
figure(13);
subplot(3,1,1)
plot(t,Pxs,'-.m',t,Pxt,'-r',t,PxToa,'-b','LineWidth',2);
legend('空间机器人X方向线动量','目标X方向线动量','复合体系统X方向线动量')
subplot(3,1,2)
plot(t,Pys,'-.m',t,Pyt,'-r',t,PyToa,'-b','LineWidth',2);
hold on;
legend('空间机器人Y方向线动量','目标Y方向线动量','复合体系统Y方向线动量')
subplot(3,1,3)
plot(t,Ls,'-.m',t,Lt,'-r',t,Ls+Lt,'-b','LineWidth',2);
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%% Lamda %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(14);
semilogx(t,Lamda,'LineWidth',2);
xlabel('Time (s)','fontsize',12);
ylabel('Variable forgetting factor','fontsize',12);
axis([0 20 -inf inf]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
time = datestr(now,29);
file = strcat('.\data\',time);
mkdir(file);
filename = strcat(file,'\data.csv');

%xlswrite('data.xls',[Bp Bv D1p D2p D3p D1v D2v D3v U1p U2p U3p U1v U2v U3v]);
% save('.\data\data.mat',[Bp Bv D1p D2p D3p D1v D2v D3v U1p U2p U3p U1v U2v U3v mT IT bT K_Error tau1_D tau2_D tau3_D tau1_U tau2_U tau3_U]);
save([file '\Torque.mat'],'tau1_D','tau2_D','tau3_D','tau1_U', 'tau2_U', 'tau3_U');  
save([file '\Error.mat'],'K_Err_out'); 
save([file '\Base_Target.mat'],'Bp', 'Bv', 'mT', 'IT','bT');
save([file '\DualArm.mat'],'D1p', 'D2p', 'D3p', 'D1v', 'D2v', 'D3v', 'U1p', 'U2p', 'U3p', 'U1v', 'U2v', 'U3v');
