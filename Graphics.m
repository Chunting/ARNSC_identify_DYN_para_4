% Parameters of the base, angles, angular velocity,
t = timer.Data;
BasePos_ = rad2deg(BasePos.Data);
BaseVel_ = rad2deg(BaseVel.Data);
BaseAcc_ = rad2deg(BaseAcc.Data);  %绕Z轴旋转的角加速度

% Angules, angular velocities of all the joints, U is right, U is lest 
JointPos_D1_ = rad2deg(JointPos_D1.Data);  % 右一关节的角度
JointPos_D2_ = rad2deg(JointPos_D2.Data);
JointPos_D3_ = rad2deg(JointPos_D3.Data);
JointVel_D1_ = rad2deg(JointVel_D1.Data);
JointVel_D2_ = rad2deg(JointVel_D2.Data);
JointVel_D3_ = rad2deg(JointVel_D3.Data);
JointPos_U1_ = rad2deg(JointPos_U1.Data);
JointPos_U2_ = rad2deg(JointPos_U2.Data);
JointPos_U3_ = rad2deg(JointPos_U3.Data);
JointVel_U1_ = rad2deg(JointVel_U1.Data);
JointVel_U2_ = rad2deg(JointVel_U2.Data);
JointVel_U3_ = rad2deg(JointVel_U3.Data);

%MassTarget_ = MassTarget.Data;
InertiaTarget_ = InertiaTarget.Data;
CoMTarget_ = CoMTarget.Data;


% Position of base in X-axis and Y-aixs
LinMomX_SR_ = LinMomX_SR.Data;
LinMomY_SR_ = LinMomY_SR.Data;
AngMom_SR_ = AngMom_SR.Data;
%AngMom_SR_(1)=0.0952;
LinMomX_Target_ = LinMomX_Target.Data;
LinMomY_Target_ = LinMomY_Target.Data;
AngMom_Target_ = AngMom_Target.Data;
%AngMom_Target_(1) = 64.9165;
LinMomX_Sum_ = LinMomX_Sum.Data;
LinMomY_Sum_ = LinMomY_Sum.Data;
AngMom_Sum_ = AngMom_Sum.Data;
n = length(AngMom_Sum_);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
JointTorque_U1_ = JointTorque_U1.Data;
JointTorque_U2_ = JointTorque_U2.Data;
JointTorque_U3_ = JointTorque_U3.Data;
JointTorque_D1_ = JointTorque_D1.Data;
JointTorque_D2_ = JointTorque_D2.Data;
JointTorque_D3_ = JointTorque_D3.Data;
C_eta = C_eta_out.Data;
Angular_Err_ = rad2deg(Angular_Err_out.Data);
K_Err = rad2deg(K_Err_out.Data);
Lamda = Lamda_out.Data;
startTime = 0.0;
stopTime = t(end);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
timeFolder = datestr(now,'yyyy_mm_dd_HH_MM');
dataPath = strcat('.\data\',timeFolder);
mkdir(dataPath);
%%%%%%%%%%%%%%%%%%%%%   1. 关节角度  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figureHandle1 = figure(1);
subplot(2,1,1)
semilogx(t,JointPos_U1_,'-b',t,JointPos_U2_,'--r',t,JointPos_U3_,'-.m','LineWidth',2);
axis([startTime stopTime -inf inf]);
hold on
% grid on
%xlabel('Time (s)','fontsize',12);
ylabel('Arm-A (deg)','fontsize',12);
hl = legend('Joint1','Joint2','Joint3','Location','northwest');
set(hl,'Orientation','horizon');
set(hl,'Box','off');
% title('Joints Angles of Arm-A');
subplot(2,1,2)
%figure(2);
semilogx(t,JointPos_D1_,'-b',t,JointPos_D2_,'--r',t,JointPos_D3_,'-.m','LineWidth',2);
axis([startTime stopTime -inf inf]);
hold on
% grid on
xlabel('Time (s)','fontsize',12);
ylabel('Arm-B (deg)','fontsize',12);
hl = legend('Joint1','Joint2','Joint3','Location','northwest');
set(hl,'Orientation','horizon');
set(hl,'Box','off');
savefig(figureHandle1,[dataPath '\JointPosition.fig']);
% title('Joints Angles of Arm-B');
%%%%%%%%%%%%%%%%%%%%%%%%%  2. 关节速度  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figureHandle3 = figure(3);
subplot(2,1,1)
semilogx(t,JointVel_U1_,'-b',t,JointVel_U2_,'--r',t,JointVel_U3_,'-.m','LineWidth',2);
% grid on
axis([startTime stopTime -inf inf]);
%xlabel('Time (s)','fontsize',12);
ylabel('Arm-A (deg/s)','fontsize',12);
hl = legend('Joint1','Joint2','Joint3','Location','northwest');
set(hl,'Orientation','horizon');
set(hl,'Box','off');
% title(' Joint rates of Arm-A');
subplot(2,1,2);
%figure(4);
semilogx(t,JointVel_D1_,'-b',t,JointVel_D2_,'--r',t,JointVel_D3_,'-.m','LineWidth',2);
axis([startTime stopTime -inf inf]);
% grid on
xlabel('Time (s)','fontsize',12);
ylabel('Arm-B (deg/s)','fontsize',12);
hl = legend('Joint1','Joint2','Joint3','Location','northwest');
set(hl,'Orientation','horizon');
set(hl,'Box','off');
savefig(figureHandle3,[dataPath '\JointVelcity.fig']);
% title(' Joint rates of Arm-B');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  3. 基座扰动  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figureHandle5 = figure(5);
subplot(2,1,1)
semilogx(t,BasePos_,'LineWidth',2);
% grid on
axis([startTime stopTime -inf inf]);
xlabel('Time (s)','fontsize',12);
ylabel('Base angles (deg)','fontsize',12);
% title(' Angles of the Base');
% BaseVel_(1)=0;
subplot(2,1,2)
semilogx(t,BaseVel_,'LineWidth',2);
% grid on
axis([startTime stopTime -inf inf]);
ylabel('Base angular velocity (deg/s)','fontsize',12);
% title(' Angular velocity of the Base');
xlabel('Time (s)','fontsize',12);
savefig(figureHandle5,[dataPath '\BaseDisturbance.fig']);
%%%%%%%%%%%%%%%%%%%%%%%%%%   4.辨识结果   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figureHandle6 = figure(6);
subplot(3,1,1)
semilogx(t,MassTarget_,'LineWidth',2);
axis([startTime stopTime 0 150]);
ylabel('Mass (kg)','fontsize',12);

subplot(3,1,2)
semilogx(t,CoMTarget_,'LineWidth',2);
ylabel('Mass center (m)','fontsize',12);
axis([startTime stopTime 0 0.5]);

subplot(3,1,3);
semilogx(t,InertiaTarget_,'LineWidth',2);
axis([startTime stopTime -inf inf]);
xlabel('Time (s)','fontsize',12);
ylabel('Inertia (kg.m^2)','fontsize',12);
savefig(figureHandle6,[dataPath '\IdentificationResult.fig']);
%%%%%%%%%%%%%%%%%%%%%%   5. 关节力矩     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figureHandle7 = figure(7);
subplot(2,1,1)
semilogx(t,JointTorque_U1_,'-b',t,JointTorque_U2_,'--r',t,JointTorque_U3_,'-.m','LineWidth',2);
axis([startTime stopTime -inf inf]);
xlabel('Time (s)','fontsize',12);
ylabel('Arm-A torque (N \cdot M)','fontsize',12);
hl = legend('\tau1','\tau2','\tau3','Location','northwest');
set(hl,'Orientation','horizon');
set(hl,'Box','off');
%savefig(figureHandle7,[dataPath '\JointTorque_U.fig']);
% title(' Joint Torque of arm-a');
subplot(2,1,2);
%figureHandle8 = figure(8);
semilogx(t,JointTorque_D1_,'-b',t,JointTorque_D2_,'--r',t,JointTorque_D3_,'-.m','LineWidth',2);
axis([startTime stopTime -inf inf]);
xlabel('Time (s)','fontsize',12);
ylabel('Arm-B torque (N \cdot M)','fontsize',12);
hl = legend('\tau1','\tau2','\tau3','Location','northwest');
set(hl,'Orientation','horizon');
set(hl,'Box','off');
savefig(figureHandle7,[dataPath '\JointTorque_D.fig']);
%%%%%%%%%%%%%%%%%%%%%    6. 动量守恒    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figureHandle9 = figure(9);
l = length(t);
plot(t,AngMom_SR_,'-.m',t,AngMom_Target_,'-r',t,AngMom_Sum_,'-b','LineWidth',2);
axis([startTime stopTime -inf inf]);
hl = legend('AngMomSR','AngMomTarget','AngMomSum');
set(hl,'Orientation','horizon');
set(hl,'Box','off');
xlabel('Time (s)','fontsize',12);
ylabel('Angular Momentum (kg·m/s)','fontsize',12);
savefig(figureHandle9,[dataPath '\AngularMomentum.fig']);
%%%%%%%%%%%%%%%%%%%%%%  7. K 误差  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figureHandle10 = figure(10);
subplot(2,1,1)
semilogx(t,K_Err(:,1:3),'LineWidth',2);
axis([startTime stopTime -inf inf]);
xlabel('Time (s)','fontsize',12);
ylabel('Arm-A','fontsize',12);
hl = legend('Joint1','Joint2','Joint3','Location','northwest');
set(hl,'Orientation','horizon');
set(hl,'Box','off');
%savefig(figureHandle10,[dataPath '\JointVel_Error_U.fig']);
%figureHandle11 = figure(11);
subplot(2,1,2)
semilogx(t,K_Err(:,4:6),'LineWidth',2);
axis([startTime stopTime -inf inf]);
xlabel('Time (s)','fontsize',12);
ylabel('Arm-B ','fontsize',12);
hl = legend('Joint1','Joint2','Joint3','Location','northwest');
set(hl,'Orientation','horizon');
set(hl,'Box','off');
savefig(figureHandle10,[dataPath '\K_Error.fig']);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figureHandle12 = figure(12);
subplot(2,1,1);
semilogx(t,Angular_Err_(:,1:3),'Linewidth',2);
axis([startTime stopTime -inf inf]);
% title('Angular Error of Joint A1 and B1');
hl = legend('A1','A2','A3','Location','northwest');
set(hl,'Orientation','horizon');
set(hl,'Box','off');
xlabel('Time (s)','fontsize',12);
ylabel('Arm-A (deg/s)','fontsize',12);
subplot(2,1,2);
semilogx(t,Angular_Err_(:,4:6),'Linewidth',2);
axis([startTime stopTime -inf inf]);
% title('Angular Error of Joint A1 and B1');
hl = legend('B1','B2','B3','Location','northwest');
set(hl,'Orientation','horizon');
set(hl,'Box','off');
xlabel('Time (s)','fontsize',12);
ylabel('Arm-B (deg/s)','fontsize',12);
savefig(figureHandle12,[dataPath '\JointAngular_Error.fig']);
%%%%%%%%%%%%%%%%%%%  9. 线动量  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
figure(13);
subplot(3,1,1)
semilogx(t,LinMomX_SR_,'-.m',t,LinMomX_Target_,'-r',t,LinMomX_Sum_,'-b','LineWidth',2);
legend('空间机器人X方向线动量','目标X方向线动量','复合体系统X方向线动量')
subplot(3,1,2)
semilogx(t,LinMomY_SR_,'-.m',t,LinMomY_Target_,'-r',t,LinMomY_Sum_,'-b','LineWidth',2);
hold on;
legend('空间机器人Y方向线动量','目标Y方向线动量','复合体系统Y方向线动量')
subplot(3,1,3)
semilogx(t,AngMom_SR_,'-.m',t,AngMom_Target_,'-r',t,AngMom_SR_+AngMom_Target_,'-b','LineWidth',2);
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%% Lamda %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figureHandle14 = figure(14);
loglog(t,Lamda,'LineWidth',2);
xlabel('Time (s)','fontsize',12);
ylabel('Variable forgetting factor','fontsize',12);
axis([startTime stopTime -inf 1]);
savefig(figureHandle14,[dataPath '\Lamda.fig']);


%xlswrite('data.xls',[BasePos_ BaseVel_ JointPos_D1_ JointPos_D2_ JointPos_D3_ JointVel_D1_ JointVel_D2_ JointVel_D3_ JointPos_U1_ JointPos_U2_ JointPos_U3_ JointVel_U1_ JointVel_U2_ JointVel_U3_]);
% save('.\data\data.mat',[BasePos_ BaseVel_ JointPos_D1_ JointPos_D2_ JointPos_D3_ JointVel_D1_ JointVel_D2_ JointVel_D3_ JointPos_U1_ JointPos_U2_ JointPos_U3_ JointVel_U1_ JointVel_U2_ JointVel_U3_ MassTarget_ InertiaTarget_ CoMTarget K_Error JointTorque_D1 JointTorque_D2_ JointTorque_D3_ JointTorque_U1_ JointTorque_U2_ JointTorque_U3_]);
save([dataPath '\Timer.mat'],'t'); 
save([dataPath '\Torque.mat'],'JointTorque_D1','JointTorque_D2_','JointTorque_D3_','JointTorque_U1_', 'JointTorque_U2_', 'JointTorque_U3_');  
save([dataPath '\K_Err_out.mat'],'K_Err_out'); 
save([dataPath '\Base_Target.mat'],'BasePos_', 'BaseVel_', 'MassTarget_', 'InertiaTarget_','CoMTarget');
save([dataPath '\Identification.mat'],'AngMom_SR_', 'AngMom_Target_', 'AngMom_Sum_');
save([dataPath '\Lamda.mat'],'Lamda');
save([dataPath '\DualArm.mat'],'JointPos_D1_', 'JointPos_D2_', 'JointPos_D3_', 'JointVel_D1_', 'JointVel_D2_', 'JointVel_D3_', 'JointPos_U1_', 'JointPos_U2_', 'JointPos_U3_', 'JointVel_U1_', 'JointVel_U2_', 'JointVel_U3_');
