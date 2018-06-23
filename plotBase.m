function plotBase(timer,B_p,B_v,B_a)
t = timer.Data;
Bp = rad2deg(B_p.Data);
Bv = rad2deg(B_v.Data);
Ba = rad2deg(B_a.Data);  %绕Z轴旋转的角加速度
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  3. 基座扰动  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure('Name','Base');
subplot(2,1,1)
plot(t,Bp,'LineWidth',2);

%axis([0 20 -0.001 0.001]);
ylabel('Base angles (deg)','fontsize',12);

Bv(1)=0;
subplot(2,1,2)
plot(t,Bv,'LineWidth',2);
%axis([0 20 -0.003 0.001]);
ylabel('Base angular velocity (deg/s)','fontsize',12);
xlabel('Time (s)','fontsize',12);
end