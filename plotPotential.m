function plotPotential(timer,D1_p,D2_p,D3_p,U1_p,U2_p,U3_p)

t = timer.Data;
Theta_Umax = zeros(size(t)) + 120;
Theta_Dmax1 = zeros(size(t)) + 60;
Theta_Dmax2 = zeros(size(t)) + 120;
D1p = rad2deg(D1_p.Data);  % 右一关节的角度
D2p = rad2deg(D2_p.Data);
D3p = rad2deg(D3_p.Data);

U1p = rad2deg(U1_p.Data);
U2p = rad2deg(U2_p.Data);
U3p = rad2deg(U3_p.Data);
pu1 = 1./(Theta_Umax.^2-U1p.^2);
pu2 = 1./(Theta_Umax.^2-U2p.^2);
pu3 = 1./(Theta_Umax.^2-U3p.^2);
pd1 = 1./(Theta_Dmax1.^2-D1p.^2);
pd2 = 1./(Theta_Dmax2.^2-D2p.^2);
pd3 = 1./(Theta_Dmax2.^2-D3p.^2);
pot = pu1+pu2+pu3+pd1+pd2+pd3;
figure('Name','Potensial Function');
plot(t,pot,'-b','LineWidth',2);
%axis([0 20 -0.01 0.01]);
hl = legend('Arm-a joint1','Arm-b joint1',3);
set(hl,'Orientation','horizon');
set(hl,'Box','off');
xlabel('Time(s)','fontsize',12);
ylabel('Potential function (deg^{-2})','fontsize',12);

end