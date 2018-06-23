function plotTorque(timer,Tau1_U,Tau2_U,Tau3_U,Tau1_D,Tau2_D,Tau3_D)
t = timer.Data;
tau1_U = Tau1_U.Data;
tau2_U = Tau2_U.Data;
tau3_U = Tau3_U.Data;
tau1_D = Tau1_D.Data;
tau2_D = Tau2_D.Data;
tau3_D = Tau3_D.Data;

figure('Name','Torque');
subplot(2,1,1)
semilogx(t,tau1_U,'-b',t,tau2_U,'--r',t,tau3_U,'-.m','LineWidth',2);
axis([0 20 -0.6 0.15]);
hl = legend('\tau1','\tau2','\tau3','Location','northwest');
set(hl,'Orientation','horizon');
set(hl,'Box','off');
%set(gca,'xticklabel',0:2:20);
ylabel('Arm-a torque (N \cdot M)','fontsize',12);

subplot(2,1,2)
semilogx(t,tau1_D,'-b',t,tau2_D,'--r',t,tau3_D,'-.m','LineWidth',2);
axis([0 20 -0.5 0.15]);
%set(gca,'xticklabel',0:2:20);
xlabel('Time (s)','fontsize',12);
ylabel('Arm-b torque (N \cdot M)','fontsize',12);
hl = legend('\tau1','\tau2','\tau3','Location','northwest');
set(hl,'Orientation','horizon');
set(hl,'Box','off');
end