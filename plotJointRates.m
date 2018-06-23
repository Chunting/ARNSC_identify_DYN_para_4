function plotJointRates(timer,D1_v,D2_v,D3_v,U1_v,U2_v,U3_v)
t = timer.Data;
% Angules, angular velocities of all the joints, U is right, U is lest
D1v = rad2deg(D1_v.Data);
D2v = rad2deg(D2_v.Data);
D3v = rad2deg(D3_v.Data);

U1v = rad2deg(U1_v.Data);
U2v = rad2deg(U2_v.Data);
U3v = rad2deg(U3_v.Data);
%%%%%%%%%%%%%%%%%%%%%%%%%  2. 关节速度  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure('Name','JointRates');
subplot(2,1,1)
plot(t,U1v,'-b',t,U2v,'--r',t,U3v,'-.m','LineWidth',2);
% grid on
%axis([0 20 -0.8 0.5]);
ylabel('Arm-a (deg/s)','fontsize',12);
hl = legend('Joint1','Joint2','Joint3','Location','northwest');
set(hl,'Orientation','horizon');
set(hl,'Box','off');
% title(' Joint rates of Arm-a');
subplot(2,1,2);
plot(t,D1v,'-b',t,D2v,'--r',t,D3v,'-.m','LineWidth',2);
% grid on
xlabel('Time (s)','fontsize',12);
ylabel('Arm-b (deg/s)','fontsize',12);
hl = legend('Joint1','Joint2','Joint3','Location','northwest');
set(hl,'Orientation','horizon');
set(hl,'Box','off');
% title(' Joint rates of Arm-b');

end