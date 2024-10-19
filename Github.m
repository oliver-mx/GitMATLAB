%% git pull
system('git pull https://github.com/oliver-mx/GitMATLAB.git');
addpath('Input_DATA','Scaled_model','Unscaled_model','Output_DATA')

%% git push
msg = input('Please enter a commit message: ', 's');
system('git status');
system('git add .');
system(['git commit -m "', msg, '"']);
system('git push https://github.com/oliver-mx/GitMATLAB.git');


%% uni PC
%system('git add .'); system('git commit -m "Neue Simulation"');system('git push https://github.com/oliver-mx/GitMATLAB.git');

%% home
%system('git status');
%system('git add .');
%system('git commit -m "fmincon sqp"');
%system('git push https://github.com/oliver-mx/GitMATLAB.git');
%system('shutdown /s /t 30');
