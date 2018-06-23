 plotJointAngle(timer,D1_p,D2_p,D3_p,U1_p,U2_p,U3_p);
 plotJointRates(timer,D1_v,D2_v,D3_v,U1_v,U2_v,U3_v);
 plotTorque(timer,Tau1_U,Tau2_U,Tau3_U,Tau1_D,Tau2_D,Tau3_D);
 plotEstimation(timer,mt,It,bt);
 plotBase(timer,B_p,B_v,B_a);
 plotMomentum(timer,L_s,L_t,L_Toa);
 plotAngularError(timer,Angular_Err_out);
 plotPotential(timer,D1_p,D2_p,D3_p,U1_p,U2_p,U3_p);

file = datestr(now,30);
file = ['.\data\' file];
mkdir(file);
save([file '\Torque.mat'],'Tau1_D','Tau2_D','Tau3_D','Tau1_U', 'Tau2_U', 'Tau3_U');  
save([file '\AngularError.mat'],'Angular_Err_out'); 
save([file '\Base_Target.mat'],'B_p', 'B_v', 'mt', 'It','bt');
save([file '\DualArm.mat'],'D1_p', 'D2_p', 'D3_p', 'D1_v', 'D2_v', 'D3_v', 'U1_p', 'U2_p', 'U3_p', 'U1_v', 'U2_v', 'U3_v');
