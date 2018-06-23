function plotLinearMom(timer)
%%%%%%%%%%%%%%%%%%%  9. Ïß¶¯Á¿  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t = timer.Data;
figure('Name','LinearMomentum');
subplot(3,1,1)
plot(t(2:4001),Pxs(2:4001),'-.m',t(2:4001),Pxt(2:4001),'-r',t(2:4001),PxToa(2:4001),'-b','LineWidth',2);
legend('Ls','Lt','Ls+Lt')
subplot(3,1,2)
plot(t(2:4001),Pys(2:4001),'-.m',t(2:4001),Pyt(2:4001),'-r',t(2:4001),PyToa(2:4001),'-b','LineWidth',2);
hold on;
legend('Ls','Lt','Ls+Lt')
subplot(3,1,3)
plot(t(2:4001),Ls(2:4001),'-.m',t(2:4001),Lt(2:4001),'-r',t(2:4001),Ls(2:4001)+Lt(2:4001),'-b','LineWidth',2);
%}

end