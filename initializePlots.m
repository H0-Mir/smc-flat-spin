function initializePlots(time_vec,states_store,input_store,s_store,desired)

%%
figure('Name','phi theta','NumberTitle','off');
hold on 
box on
plot(time_vec,states_store(:,7)*180/pi,'k--','LineWidth',1.5)
plot(time_vec,states_store(:,8)*180/pi,'b-.','LineWidth',1.5)
plot(time_vec,desired(1)*180/pi*ones(length(time_vec),1),'k--','LineWidth',0.5)
plot(time_vec,desired(2)*180/pi*ones(length(time_vec),1),'b-.','LineWidth',0.5)
xlim([0 time_vec(end)])
xlabel('Time (s)')
ylabel('Angle (deg)')
legend('$\phi$','$\theta$','$\phi_d$','$\theta_d$','Interpreter','latex','NumColumns',2);
%%
psi_dot = sin(states_store(:,7))./cos(states_store(:,8)).*states_store(:,5)...
         +cos(states_store(:,7))./cos(states_store(:,8)).*states_store(:,6);
figure('Name','psi_dot','NumberTitle','off');
hold on 
box on
plot(time_vec,psi_dot*180/pi,'r:','LineWidth',1.5)
plot(time_vec,desired(3)*180/pi*ones(length(time_vec),1),'r-.','LineWidth',0.5)

xlim([0 time_vec(end)])
xlabel('Time (s)')
ylabel('Yaw rate (deg/s)')
legend('$\dot \psi$','$\dot \psi_d$','Interpreter','latex','NumColumns',2);
%%
figure('Name','Altitude','NumberTitle','off');
hold on 
box on
plot(time_vec,states_store(:,12),'k-.','LineWidth',1.5)
plot(time_vec,desired(4)*ones(length(time_vec),1),'k-','LineWidth',0.5)

xlim([0 time_vec(end)])
ylim([-100 2000])
xlabel('Time (s)')
ylabel('Altitude (m)')
legend('$h$','$h_d$','Interpreter','latex');
%%
figure('Name','Angular rate','NumberTitle','off');
hold on 
box on
plot(time_vec,states_store(:,4)*180/pi,'k--','LineWidth',1.5)
plot(time_vec,states_store(:,5)*180/pi,'b-.','LineWidth',1.5)
plot(time_vec,states_store(:,6)*180/pi,'r:','LineWidth',1.5)
xlim([0 time_vec(end)])
xlabel('Time (s)')
ylabel('Angular rate (deg/s)')
legend('$p$','$q$','$r$','Interpreter','latex','NumColumns',3);
%%
V = sqrt(sum(states_store(:,1:3).^2,2));
alpha = atan(states_store(:,3)./states_store(:,1));
beta  = asin(states_store(:,2)./V);
gamma = states_store(:,8)-alpha;

figure('Name','alpha beta gamma','NumberTitle','off');
hold on 
box on
plot(time_vec,alpha*180/pi,'k--','LineWidth',1.5)
plot(time_vec,beta*180/pi,'b-.','LineWidth',1.5)
plot(time_vec,gamma*180/pi,'r:','LineWidth',1.5)
xlim([0 time_vec(end)])
xlabel('Time (s)')
ylabel('Angle (deg)')
legend('$\alpha$','$\beta$', '$\gamma$','Interpreter','latex','NumColumns',3);
%%
figure('Name','Velocity','NumberTitle','off');
hold on 
box on
plot(time_vec,V,'LineWidth',1.5,'Color',[1 0 1])
plot(time_vec,states_store(:,1),'k--','LineWidth',1.5)
plot(time_vec,states_store(:,2),'b-.','LineWidth',1.5)
plot(time_vec,states_store(:,3),'r:','LineWidth',1.5)
xlim([0 time_vec(end)])
xlabel('Time (s)')
ylabel('Velocity (m/s)')
legend('$V$','$u$', '$v$' ,'$w$','Interpreter','latex','NumColumns',4);
%% %%%%%%%%%%%%%%%%%%%%

figure('Name','Elevator','NumberTitle','off');
plot(time_vec,input_store(:,1)*180/pi,'k--','LineWidth',1.5)
xlim([0 time_vec(end)])
ylabel('Elevator input (deg)')
xlabel('Time (s)')

figure('Name','Aileron','NumberTitle','off');
plot(time_vec,input_store(:,2)*180/pi,'b-.','LineWidth',1.5)
xlim([0 time_vec(end)])
ylabel('Aileron input (deg)')
xlabel('Time (s)')

figure('Name','Rudder','NumberTitle','off');
plot(time_vec,input_store(:,3)*180/pi,'r:','LineWidth',1.5)
xlim([0 time_vec(end)])
ylabel('Rudder input (deg)')
xlabel('Time (s)')

figure('Name','Throttle','NumberTitle','off');
plot(time_vec,input_store(:,4),'b','LineWidth',1.5)
xlim([0 time_vec(end)])
ylabel('Throttle')
xlabel('Time (s)')

%%
figure('Name','Sliding surface 1,2,3','NumberTitle','off');
hold on
box on
plot(time_vec(1:end-1),s_store(:,1)*180/pi,'r:','LineWidth',1.5)
plot(time_vec(1:end-1),s_store(:,2)*180/pi,'b-.','LineWidth',1.5)
plot(time_vec(1:end-1),s_store(:,3)*180/pi,'k--','LineWidth',1.5)
xlim([0 time_vec(end)])
ylabel('Sliding surface')
xlabel('Time (s)')
legend('$s_1$','$s_2$', '$s_3$','Interpreter','latex','NumColumns',3);
%%
figure('Name','Sliding surface 4','NumberTitle','off');
hold on
box on
plot(time_vec(1:end-1),s_store(:,4),'r:','LineWidth',1.5)
xlim([0 time_vec(end)])
ylabel('Sliding surface')
xlabel('Time (s)')
legend('$s_4$','Interpreter','latex');
%%
figure('Name','trajectory','NumberTitle','off');
hold on
plot3(states_store(:,10),states_store(:,11),states_store(:,12),'b')
xlabel('X (m)');ylabel('Y (m)'); zlabel('h (m)');
box on
view(70,30)
axis equal
