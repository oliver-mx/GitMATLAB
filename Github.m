%% git pull
system('git pull https://github.com/oliver-mx/GitMATLAB.git');
addpath('Input_DATA','Scaled_model','Unscaled_model','Output_DATA')

%% git push
msg = input('Please enter a commit message: ', 's');
system('git status');
system('git add .');
system(['git commit -m "', msg, '"']);
system('git push https://github.com/oliver-mx/GitMATLAB.git');

%system('git add .'); system('git commit -m "Neue Simulationen"');system('git push https://github.com/oliver-mx/GitMATLAB.git');



