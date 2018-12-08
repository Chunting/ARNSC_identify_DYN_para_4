
function Y = Runge_Kutta4(t,ytemp,dfunc,h,u)
% disp "================================================="
% ytemp;
K1=dfunc(t,ytemp,u);
K2=dfunc(t+h/2,ytemp+h/2*K1,u);
K3=dfunc(t+h/2,ytemp+h/2*K2,u);
K4=dfunc(t+h,ytemp+h*K3,u);
Y=ytemp+h/6*(K1+2*K2+2*K3+K4);
% disp "================================================="
end