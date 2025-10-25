clc;  close all;

%% Parameters
load("params.mat")

desired = [0; deg2rad(17.12); 0; 0;
           0; 0; 0; 0];
init_conds = [      [20.6;  3.1;  56.4];   
            deg2rad([-22.92; 2.86; -83.1]); 
            deg2rad([-1.15; -19.25; 0]); 
                    [0;  0;  1557]          ];
input0 = [deg2rad(-24.1); 0;  0; 0.4];
states_store = init_conds';
time_vec = 0;
input_store = input0';
s_store = [];
t_f = 100;
dt  = 0.1;
%% Simulation
fprintf('Starting Simulation...\n')
% Progress bar
dots = '....................';
bars = '||||||||||||||||||||';
rewind = '\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b';
fprintf(strcat('Simulating: ',dots));
progress_mark = 0;
progress_epoch = 0;


tic
for t = 0:dt:t_f-dt
    [ctrl_input,s_out]  = Control(params,ctrl_params,states_store(end, :)',input_store(end,:),desired);
    if  t<10.5
         input = input0;
    else
         input = ctrl_input;
    end
    [~, states_temp] = ode45(@(t, states) Aircraft_EOM(t,params,states,input), [t t+dt], states_store(end, :));
    states_store = [states_store; states_temp(end, :)]; 
    time_vec = [time_vec; t + dt];          
    input_store = [input_store; input'];
    s_store = [s_store; s_out'];
        % Update progress bar
    if (t - progress_epoch) > (t_f/20)
        progress_mark = progress_mark + 1;
        progress_epoch = t;
        fprintf(strcat(rewind,bars(1:progress_mark),...
            dots(1:(20 - progress_mark))));
    end 
    
end
% Complete progress bar
fprintf(strcat(rewind,bars,'\n'));
fprintf('Simulation ends.\n')
toc
%% Plotting
initializePlots(time_vec,states_store,input_store,s_store,desired)