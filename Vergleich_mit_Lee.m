%% Basic simulations for all cases (ideal, ICP, ICP+ECP)
% press [ctrl]+[enter] to run code sections
addpath('Input_DATA','Scaled_model','Unscaled_model','Output_DATA')

%% quick test
clc;
%[a3,b3]=fun_scaled([-0.335; 0 ; 10 ; 1.000009; 0.063],0,'fig',1e4,1e-6);
[a3,b3]=fun_scaled([-0.335; 0 ; 5 ; 1.05; 0.063],0,'sol',1e4,1e-6);

ev(a3,[5 6])








%% Vergleich mit Senthil Table 2
clc
output=fun_unscaled([0; 50.47; 48.23],0,'fig',1e4,1e-3);
disp('Operating cond:  P_d(0) = 50.47 bar      P_d(L) = 48.23 bar')
disp(' ')
disp('Experiment:      Recovery = 33.1 %       C_Permeate = 44 ppm')
disp(['Our model:       Recovery = ',num2str(output(4)),' %    C_Permeate = ',num2str(output(9)),' ppm   ==>   J_d(0) = ',num2str(131.7146216211166*output(6)),' kg/sm'])
str1=100*(output(4)-33.1)/33.1;str2=100*(output(9)-44)/44;
if str1 > 0; str1=['+',num2str(str1)]; else str1=num2str(str1); end
if str2 > 0; str2=['+',num2str(str2)]; else str2=num2str(str2); end
warning(['                   ',str1,' %                ',str2,' %'])
disp('_______________')
J_d_in=131.7146216211166*output(6);
output=fun_unscaled([J_d_in; 55.81; 0],0,'sol',1e4,1e-3);
disp(['Operating cond:  P_d(0) = 55.81 bar      J_d(0) = ',num2str(131.7146216211166*output(6)),' kg/sm'] )
disp(' ')
disp('Experiment:      Recovery = 39.6 %       C_Permeate = 49 ppm         P_d(L) = 54.72')
disp(['Our model:       Recovery = ',num2str(output(4)),' %     C_Permeate = ',num2str(output(9)),' ppm    P_d(L) = ',num2str(output(11)),' bar'])
str1=100*(output(4)-39.6)/39.6;str2=100*(output(9)-49)/49;str3=100*(output(11)-54.72)/54.72;
if str1 > 0; str1=['+',num2str(str1)]; else str1=num2str(str1); end
if str2 > 0; str2=['+',num2str(str2)]; else str2=num2str(str2); end
if str3 > 0; str3=['+',num2str(str3)]; else str3=num2str(str3); end
warning(['                   ',str1,' %                 ',str2,' %               ',str3,' %'])
disp('_______________')
output=fun_unscaled([J_d_in; 60.28; 0],0,'sol',1e4,1e-3);
disp(['Operating cond:  P_d(0) = 60.28 bar      J_d(0) = ',num2str(131.7146216211166*output(6)),' kg/sm'] )
disp(' ')
disp('Experiment:      Recovery = 44.5 %       C_Permeate = 52 ppm         P_d(L) = 59.22')
disp(['Our model:       Recovery = ',num2str(output(4)),' %    C_Permeate = ',num2str(output(9)),' ppm    P_d(L) = ',num2str(output(11)),' bar'])
str1=100*(output(4)-44.5)/44.5;str2=100*(output(9)-52)/52;str3=100*(output(11)-59.22)/59.22;
if str1 > 0; str1=['+',num2str(str1)]; else str1=num2str(str1); end
if str2 > 0; str2=['+',num2str(str2)]; else str2=num2str(str2); end
if str3 > 0; str3=['+',num2str(str3)]; else str3=num2str(str3); end
warning(['                   ',str1,' %                 ',str2,' %               ',str3,' %'])
disp('_______________')

