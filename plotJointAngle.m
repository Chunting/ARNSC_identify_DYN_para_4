function plotJointAngle(timer,D1_p,D2_p,D3_p,U1_p,U2_p,U3_p)

t = timer.Data;
a=0.6;
Theta_Umax = zeros(size(t)) + 120;
Theta_Dmax1 = zeros(size(t)) + 60;
aTheta_Dmax1 = a*Theta_Dmax1;
aTheta_Umax = a*Theta_Umax;
% Theta_Dmax2 = zeros(size(t)) + 120;
% Angules, angular velocities of all the joints, U is right, U is lest
D1p = rad2deg(D1_p.Data);  % 右一关节的角度
D2p = rad2deg(D2_p.Data);
D3p = rad2deg(D3_p.Data);
U1p = rad2deg(U1_p.Data);
U2p = rad2deg(U2_p.Data);
U3p = rad2deg(U3_p.Data);

%%%%%%%%%%%%%%%%%%%%%   1. 关节角度  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure('Name','JointAngles_1');
subplot(2,1,1)
plot(t,U1p,'-b',t,U2p,'--r',t,U3p,'-.m','LineWidth',2);
hold on

%plot(t,Theta_Umax,':k',t,Theta_Umin,':k','LineWidth',0.5);

% grid on
ylabel('Arm-a (deg)','fontsize',12);
hl = legend('Joint1','Joint2','Joint3','Location','northwest');
set(hl,'Orientation','horizon');
set(hl,'Box','off');

subplot(2,1,2)
plot(t,D1p,'-b',t,D2p,'--r',t,D3p,'-.m','LineWidth',2);
hold on
%plot(t,Theta_Dmax1,':k',t,Theta_Dmin1,':k','LineWidth',0.5);

% axis([0 20 -2 2]);
xlabel('Time (s)','fontsize',12);
ylabel('Arm-b (deg)','fontsize',12);
hl = legend('Joint1','Joint2','Joint3','Location','northwest');
set(hl,'Orientation','horizon');
set(hl,'Box','off');

figure('Name','JointAngles_2');
plot(t,D2p,'--r',t,D3p,'-.m','LineWidth',2);
hold on;
plot(t,2*Theta_Dmax1,'-k',t,2*aTheta_Dmax1,'--g',t,-2*Theta_Dmax1,'-k',t,-2*aTheta_Dmax1,'--g','LineWidth',2);

axis([0 20 -150 150]);
xlabel('Time (s)','fontsize',12);
ylabel('Arm-b joint angles (deg)','fontsize',12);
hl = legend('joint2','joint3','Location','northwest');
set(hl,'Box','off','Orientation','horizon');

figure('Name','JointAngles_3');
plot(t,U1p,'-b',t,U2p,'--r',t,U3p,'-.m','LineWidth',2);
hold on
plot(t,Theta_Umax,'-k',t,aTheta_Umax,'--g',t,-Theta_Umax,'-k',t,-aTheta_Umax,'--g','LineWidth',2);
%hl = legend('Arm-b joint1','Limit','Safe area','Location','northwest');
%axis([0 20 -100 100]);
xlabel('Time (s)','fontsize',12);
ylabel('Arm-a joint angles (deg)','fontsize',12);
set(hl,'Box','off','Orientation','horizon');
end