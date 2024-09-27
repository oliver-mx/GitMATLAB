%%
addpath('Input_DATA','Scaled_model','Unscaled_model','Output_DATA')
%% test with ideal membrane
clc,
[a,b,c]=fun_test([55e5;0],0,'fig',1e5,1e-3);
display(' ')
display('Test (without ERD):')
display(a(1:4))
display(' ')
display('Senthil (with ERD):')
display('   -2.2000  419.0400      ?      38.0000')
%
% 0.291 feed flow rate
% 55bar inlet pressure
%
% 46g/m^3 product water concentration
% product water flow rate > 0.1164 m^3/s
%

%% unscaeld test:
clc,

%[a,b,c]=fun_unscaled([55.81e5;0e5],0,'fig',1e4,1e-3);
[a,b,c]=fun_scaled([55.81;54.72],0,'fig',1e4,1e-3);
disp('Scaled:')
disp(['   ',num2str(a(1)),'   ',num2str(a(2))])
disp(c)
[a,b,c]=fun_unscaled([55.81e5;54.72e5],0,'fig',1e4,1e-3);
disp('Unscaled:')
disp(['   ',num2str(a(1)),'   ',num2str(a(2))])
disp(c)
disp('Senthil (with ERD):')
disp('   -2.2000  419.0400')
disp('    0.2910    0.1160   39.6000   49.0000')

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






