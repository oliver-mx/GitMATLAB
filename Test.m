%% press [ctrl]+[enter] to run code sections
addpath('Input_DATA','Scaled_model','Unscaled_model','Output_DATA')

%% unscaeld test:
clc,
%
[a1,b1,c1]=fun_scaled(  [55.81e5;54.72e5],.1,'fig',1e4,1e-3);
[a2,b2,c2]=fun_unscaled([55.81e5;54.72e5],.1,'fig',1e4,1e-3);
%
disp('%             Scaled and unscaled')
disp(['% SEC_net     ',num2str(round(1e4*a1(1))/1e4),'    ',num2str(round(1e4*a2(1))/1e4)])
disp(['% FW           ',num2str(round(1e4*a1(2))/1e4),'     ',num2str(round(1e4*a2(2))/1e4)])
disp(['% REC         ',num2str(round(1e4*a1(4))/1e4),'    ',num2str(round(1e4*a2(4))/1e4)])
disp(['% C (in ppm)  ',num2str(c1(4)),'    ',num2str(c2(4))])

%% ideal
%             Scaled and unscaled
% SEC_net     -3.4004    -3.4
% FW           0.8095     0.7982
% REC         48.9183    48.9253
% C (in ppm)   0          0

%% nonideal (fixed beta)
%             Scaled and unscaled
% SEC_net     -3.5037    -3.5037
% FW           0.6716     0.6716 
% REC         47.3632    47.3629
% C (in ppm)  51.9032    51.9027

%% ICP+ECP (fixed beta)
% SEC_net     -3.5228    -3.5228
% FW           0.6615     0.6615
% REC         47.0657    47.0654
% C (in ppm)  70.5322    70.5317

%% testing unscaled SWRO with ERD
clc
[a3,b3,c3]=fun_unscaled([55.81e5;54.72e5],.2,'fig',1e4,1e-3);
disp(' ')
disp('simple RO')
disp(a2)
disp(b2)
disp('SWRO+ERD')
disp(a3)
disp(b3)















%%
a2=[NaN, NaN, NaN, NaN];
c2=[NaN, NaN, NaN, NaN];
%disp('Senthil (with ERD):')
%disp('   -2.2000  419.0400')
%disp('    0.2910    0.1160   39.6000   49.0000')

% 0.291 feed flow rate
% 55bar inlet pressure
%
% 46g/m^3 product water concentration
% product water flow rate > 0.1164 m^3/s
%

%% Senthil experimental data
clc,
% test 1 
[a,b,c]=fun_unscaled([50.47e5;48.23e5],0,'fig',1e4,1e-3);
disp('----------------------------------------')
disp('Our model:')
disp(c)
disp('Experimental data (Senthil):')
disp('    0.2910    0.0934   33.1000   44.0000')
disp(' ')
disp('----------------------------------------')
% test 2 
[a,b,c]=fun_unscaled([55.81e5;54.72e5],0,'fig',1e4,1e-3);
disp('Our model:')
disp(c)
disp('Experimental data (Senthil):')
disp('    0.2910    0.1160   39.6000   49.0000')
disp(' ')
disp('----------------------------------------')
% test 3 
[a,b,c]=fun_unscaled([60.28e5;59.22e5],0,'fig',1e4,1e-3);
disp('Our model:')
disp(c)
disp('Experimental data (Senthil):')
disp('    0.2911    0.1295   44.5000   52.0000')



%%
system('git pull https://github.com/oliver-mx/GitMATLAB.git');

%%
system('git status');
system('git add .');
system('git commit -m "Further work on ERD..."');
system('git push https://github.com/oliver-mx/GitMATLAB.git');
