%%%%%%%%%%%%%%%%%%%%%%  7. 关节角速度误差  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plotAngularError(timer,Angular_Err_out)
t = timer.Data;
Angular_Err = rad2deg(Angular_Err_out.Data);
figure('Name','AngularError_1');
subplot(2,1,1)
plot(t,Angular_Err(:,1),t,Angular_Err(:,2),t,Angular_Err(:,3),'LineWidth',2);
%axis([0 10 -1.5e-4 1.5e-4]);
set(gca,'xticklabel',0:2:20);
%set(gca,'ytickLabel',{'0×10^-3','0.5×10^-3','1×10^-3','1.5×10^-3','2×10^-3'})
ylabel('Arm-a (deg/s)','fontsize',12);
hl = legend('joint1','joint2','joint3','Location','northwest');
set(hl,'Orientation','horizon');
set(hl,'Box','off');

subplot(2,1,2);
%plot(t,Angular_Err(:,4),t,Angular_Err(:,5),t,Angular_Err(:,6),'LineWidth',2);
plot(t,Angular_Err(:,4),t,Angular_Err(:,5),t,Angular_Err(:,6),'LineWidth',2);
%axis([0 10 -1.5e-4 1.5e-4]);
set(gca,'xticklabel',0:2:20);
xlabel('Time (s)','fontsize',12);
ylabel('Arm-b (deg/s)','fontsize',12);
hl = legend('joint1','joint2','joint3','Location','northwest');
set(hl,'Orientation','horizon');
set(hl,'Box','off');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure('Name','AngularError_2');
plot(t,Angular_Err,'LineWidth',2);
hl = legend('A1','A2','A3','B1','B2','B3','Location','northwest');
set(hl,'Orientation','horizon');
set(hl,'Box','off');
xlabel('Time (s)','fontsize',12);
ylabel('Joint angular error (deg)','fontsize',12);
end