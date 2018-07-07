filePath = '.\data\Lamda\';
filePath = '.\';
K_Err_09_File = [filePath 'K_Err_out_09.mat'];
K_Err_099_File = [filePath 'K_Err_out_099.mat'];
K_Err_0999_File = [filePath 'K_Err_out_0999.mat'];
K_Err_File = [filePath 'K_Err_out.mat'];
time_File = [filePath 'Timer.mat'];
Base_Target_09_File = [filePath 'Base_Target_09.mat'];
Base_Target_099_File = [filePath 'Base_Target_099.mat'];
Base_Target_0999_File = [filePath 'Base_Target_0999.mat'];
Base_Target_File = [filePath 'Base_Target.mat'];
time_File = [filePath 'Timer.mat'];
load(time_File);

K_Err_09_Struct = load(K_Err_09_File, 'K_Err_out');
K_Err_09_array = K_Err_09_Struct.K_Err_out;
K_Err_09_data = K_Err_09_array.Data;

K_Err_099_Struct = load(K_Err_099_File, 'K_Err_out');
K_Err_099_array = K_Err_099_Struct.K_Err_out;
K_Err_099_data = K_Err_099_array.Data;

K_Err_0999_Struct = load(K_Err_0999_File, 'K_Err_out');
K_Err_0999_array = K_Err_0999_Struct.K_Err_out;
K_Err_0999_data = K_Err_0999_array.Data;

K_Err_Struct = load(K_Err_File, 'K_Err_out');
K_Err_array = K_Err_Struct.K_Err_out;
K_Err_data = K_Err_array.Data;

startTime = 0.0;
stopTime = t(end);

varBasePos = 'BasePos_';
varBaseVel = 'BaseVel_';
Base_Target_09_Struct = load(Base_Target_09_File, varBasePos, varBaseVel);
BasePos_09_data = Base_Target_09_Struct.(varBasePos);
BaseVel_09_data = Base_Target_09_Struct.(varBaseVel);

Base_Target_099_Struct = load(Base_Target_099_File, varBasePos, varBaseVel);
BasePos_099_data = Base_Target_099_Struct.(varBasePos);
BaseVel_099_data = Base_Target_099_Struct.(varBaseVel);

Base_Target_0999_Struct = load(Base_Target_0999_File, varBasePos, varBaseVel);
BasePos_0999_data = Base_Target_0999_Struct.(varBasePos);
BaseVel_0999_data = Base_Target_0999_Struct.(varBaseVel);

Base_Target_Struct = load(Base_Target_File, varBasePos, varBaseVel);
BasePos_data = Base_Target_Struct.(varBasePos);
BaseVel_data = Base_Target_Struct.(varBaseVel);

figureHandleBaseVel = figure(1);
semilogx(t,BaseVel_09_data,'b', ...
         t,BaseVel_099_data,'r',...
         t,BaseVel_0999_data,'m',...
         t,BaseVel_data,'g','LineWidth',2);
axis([startTime stopTime -0.01 inf]);
ylabel('Base angular velocity (deg/s)','fontsize',12);
xlabel('Time (s)','fontsize',12);
hl = legend('\lambda=0.9','\lambda=0.99','\lambda=0.999','VFF','Location','northwest');
savefig(figureHandleBaseVel,[filePath '\BaseDisturbance_Vel.fig']);

figureHandleBasePos = figure(2);
semilogx(t,BasePos_09_data,'b', ...
         t,BasePos_099_data,'r',...
         t,BasePos_0999_data,'m',...
         t,BasePos_data,'g','LineWidth',2);
axis([startTime stopTime -0.05 inf]);
ylabel('Base angles (deg)','fontsize',12);
xlabel('Time (s)','fontsize',12);
hl = legend('\lambda=0.9','\lambda=0.99','\lambda=0.999','VFF','Location','northwest');
savefig(figureHandleBasePos,[filePath '\BaseDisturbance_Pos.fig']);
%{
figureHandle10 = figure(10);
subplot(2,1,1)
semilogx(t,K_Err_09_data(:,1:3),'LineWidth',2);
axis([startTime stopTime -inf inf]);
xlabel('Time (s)','fontsize',12);
ylabel('Arm-A','fontsize',12);
hl = legend('Joint1','Joint2','Joint3','Location','northwest');
set(hl,'Orientation','horizon');
set(hl,'Box','off');
%savefig(figureHandle10,[dataPath '\JointVel_Error_U.fig']);
%figureHandle11 = figure(11);
subplot(2,1,2)
semilogx(t,K_Err_09_data(:,4:6),'LineWidth',2);
axis([startTime stopTime -inf inf]);
xlabel('Time (s)','fontsize',12);
ylabel('Arm-B ','fontsize',12);
hl = legend('Joint1','Joint2','Joint3','Location','northwest');
set(hl,'Orientation','horizon');
set(hl,'Box','off');


for i = 1:6
    K_Err_figureHandle = figure(i);
    semilogx(t,K_Err_09_data(:,i),'b', t,K_Err_099_data(:,i),'r', t,K_Err_0999_data(:,i),'m',t,K_Err_data(:,i),'g','Linewidth',2);
    axis([startTime stopTime -inf inf]);
    xlabel('Time (s)','fontsize',12);
    ylabel(['Arm-A Joint ' int2str(i)],'fontsize',12);
    hl = legend('\lambda=0.9','\lambda=0.99','\lambda=0.999','VFF','Location','northwest');
    set(hl,'Orientation','horizon');
    set(hl,'Box','off');
    figName = [filePath 'K_Err_Joint_' int2str(i) '.fig'];
    savefig(K_Err_figureHandle,figName);
end
%}
