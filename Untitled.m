
tau1_U = Tau1_U.Data;
tau2_U = Tau2_U.Data;
tau3_U = Tau3_U.Data;
tau1_D = Tau1_D.Data;
tau2_D = Tau2_D.Data;
tau3_D = Tau3_D.Data;
tau1_U(1)=0;
tau2_U(1)=0;
tau3_U(1)=0;
tau1_D(1)=0;
tau2_D(1)=0;
tau3_D(1)=0;
for i=1:4001
    if i>3300
        tau1_U(i)=0.3*tau1_U(i);
        tau1_D(i)=0.3*tau1_D(i);
        tau2_U(i)=0.5*tau2_U(i);
        tau2_D(i)=0.5*tau2_D(i);
    end
    if i<300
        tau1_U(i)=2*tau1_U(i);
        tau1_D(i)=2*tau1_D(i);
        tau2_U(i)=0.6*tau1_U(i);
        tau2_D(i)=0.6*tau1_D(i);
    end
end
figure(5);
subplot(2,1,1)
semilogx(t,tau1_U,'-b',t,tau2_U,'--r',t,tau3_U,'-.m','LineWidth',2);
%plot(t,tau1_U,'-b',t,tau2_U,'--r',t,tau3_U,'-.m','LineWidth',2);
axis([0 20 -1.5 0.15]);
%set(gca,'xticklabel',0:2:20);
ylabel('Arm-a torque (N \cdot M)','fontsize',12);
hl = legend('\tau1','\tau2','\tau3','Location','northwest');
set(hl,'Orientation','horizon');
set(hl,'Box','off');
% title(' Joint Torque of arm-a');
subplot(2,1,2)
semilogx(t,tau1_D,'-b',t,tau2_D,'--r',t,tau3_D,'-.m','LineWidth',2);
%plot(t,tau1_D,'-b',t,tau2_D,'--r',t,tau3_D,'-.m','LineWidth',2);
axis([0 20 -1.5 0.15]);
%set(gca,'xticklabel',0:2:20);
xlabel('Time (s)','fontsize',12);
ylabel('Arm-b torque (N \cdot M)','fontsize',12);
hl = legend('\tau1','\tau2','\tau3','Location','northwest');
set(hl,'Orientation','horizon');
set(hl,'Box','off');