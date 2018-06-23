%%%%%%%%%%%%%%%%%%%%%    6. ∂Ø¡ø ÿ∫„    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plotMomentum(timer,L_s,L_t,L_Toa)
t = timer.Data;
Ls = L_s.Data;
Lt = L_t.Data;
LToa = L_Toa.Data;
n = length(LToa);

figure('Name','AngularMomentum');
plot(t(2:4001),Ls(2:4001),'-.m',t(2:4001),Lt(2:4001),'-r',t(2:4001),LToa(2:4001),'-b','LineWidth',2);
% % grid on
% axis([0 20 -5 95]);
% title('Angular Momentum');
hl = legend('L_m  ','L_t  ','L_m+L_t',1);
set(hl,'Orientation','horizon');
set(hl,'Box','off');
xlabel('Time (s)','fontsize',12);
ylabel('Angular Momentum (kg°§m/s)','fontsize',12);

end