function plotEstimation(timer,mt,It,bt)
t = timer.Data;
mT = mt.Data;
IT = It.Data;
bT = bt.Data;
figure('Name','Estimation');

subplot(3,1,1)
semilogx(t,mT,'LineWidth',2)
axis([0 20 0 300]);
ylabel('Mass (kg)','fontsize',12);


subplot(3,1,2)
semilogx(t,bT,'LineWidth',2);
ylabel('Mass center (m)','fontsize',12);
axis([0 20 -0.5 0.5]);

subplot(3,1,3);
semilogx(t,IT,'LineWidth',2);
axis([0 20 -100 300]);
xlabel('Time (s)','fontsize',12);
ylabel('Inertia (kg.m^2)','fontsize',12);
end
