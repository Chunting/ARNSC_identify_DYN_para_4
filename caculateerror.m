t = timer.Data;
Angular_Err = rad2deg(Angular_Err_out.Data);
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
etU1(1) = 0;
etU2(1) = 0;
etU3(1) = 0;
etD1(1) = 0;
etD2(1) = 0;
etD3(1) = 0;
%{
for i=2:2001
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
    etU1(i) = (tau1_U(i)-tau1_U(i-1))*0.01;
    etU2(i) = (tau2_U(i)-tau2_U(i-1))*0.01;
    etU3(i) = (tau3_U(i)-tau3_U(i-1))*0.01;
    etD1(i) = (tau1_D(i)-tau1_D(i-1))*0.01;
    etD2(i) = (tau2_D(i)-tau2_D(i-1))*0.01;
    etD3(i) = (tau3_D(i)-tau3_D(i-1))*0.01;
    
end

tb=6;
figure('Name','Error');
subplot(2,1,1)
semilogx(t,etU1,'-b',t,etU2,'--r',t,etU3,'-.m','LineWidth',2);
%plot(t,tau1_U,'-b',t,tau2_U,'--r',t,tau3_U,'-.m','LineWidth',2);
axis([0 20 -1e-2 1e-3]);
%set(gca,'xticklabel',0:2:20);
ylabel('Arm-a torque (N \cdot M)','fontsize',12);
hl = legend('\tau1','\tau2','\tau3','Location','northwest');
set(hl,'Orientation','horizon');
set(hl,'Box','off');
% title(' Joint Torque of arm-a');
subplot(2,1,2)
semilogx(t,etD1,'-b',t,etD2,'--r',t,etD3,'-.m','LineWidth',2);
%plot(t,tau1_D,'-b',t,tau2_D,'--r',t,tau3_D,'-.m','LineWidth',2);
axis([0 20 -1e-3 1e-3]);
%set(gca,'xticklabel',0:2:20);
xlabel('Time (s)','fontsize',12);
ylabel('Arm-b torque (N \cdot M)','fontsize',12);
hl = legend('\tau1','\tau2','\tau3','Location','northwest');
set(hl,'Orientation','horizon');
set(hl,'Box','off');

a = 100;
%}
figure('Name','AngularError');
subplot(2,1,1)

semilogx(t,Angular_Err(:,1),t,Angular_Err(:,2),t,Angular_Err(:,3),'LineWidth',2);
axis([0 20 -1.5 0.6]);
ylabel('Arm-a (deg/s)','fontsize',12);
hl = legend('joint1','joint2','joint3','Location','northwest');
set(hl,'Orientation','horizon');
set(hl,'Box','off');


subplot(2,1,2);
semilogx(t,Angular_Err(:,4),t,Angular_Err(:,5),t,Angular_Err(:,6),'LineWidth',2);
axis([0 20 -1 0.5]);
%set(gca,'xticklabel',0:2:20);
xlabel('Time (s)','fontsize',12);
ylabel('Arm-b (deg/s)','fontsize',12);
hl = legend('joint1','joint2','joint3','Location','northwest');
set(hl,'Orientation','horizon');
set(hl,'Box','off');
