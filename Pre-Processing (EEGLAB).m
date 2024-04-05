%% This scripts can be used via EEGLAB Interface. Download EEGLAB: https://sccn.ucsd.edu/eeglab/index.php

%% Pre processing
% Line#1 for loading .BDF file
% Line#2 for loading .set file
% Line#3 for filtering data 
% Line#4 for changing the sampling rate



% Line#1 Specify the folder containing the file, ex: C:\---
EEG = pop_biosig('C:\ * FOLDER LOCATION * .bdf');


% Line#2 Specify filename and folder containing the file
EEG = pop_loadset('filename','Subject_17_ICA.set','filepath','C:\ * FOLDER LOCATION *');


% Line#3 Specify low cutoff and High cutoff of the pass band
EEG = pop_eegfiltnew(EEG, 'locutoff',0.1,'hicutoff',30,'plotfreqz',1);
     

% Line#4 Downsampling 
EEG = pop_resample( EEG, 128);


% Line#5 Running ICA 
EEG = pop_runica(EEG, 'icatype', 'runica', 'extended',1,'interrupt','on');

% Line#6 Channel location 
EEG=pop_chanedit(EEG ,'load',{'C:\ * FOLDER LOCATION * .elp' 'filetype' 'autodetect'});
  
% Line#7 Display components with button to visualize their properities after ICA
pop_selectcomps(EEG, [1:30] );
 
% Line#8 Running MARA for artifact rejection, gives probability of channel being an artifact      
pop_processMARA ( ALLEEG,EEG,CURRENTSET )
        
% Line#9 Removing channels (bad)
EEG = pop_subcomp( EEG, [1    6   18   30   31   33   55   59   61   73   74   81   83   86   87  104  109  117  119  122], 0);
     
% Line#10 Extracting epochs (EEG, events, timelimits)
EEG = pop_epoch( EEG, {  }, [-0.2         1.8], 'newname', 'ICA epochs', 'epochinfo', 'yes');
     
% Line#11 Removing baseline
EEG = pop_rmbase(EEG, [-200 0],0);

% Line#12 Rejecting trials
pop_eegplot( EEG, 1, 1, 1);

% Line#13 Rejecting epoch
EEG = pop_rejepoch( EEG, [1 2] ,0);

% Line#14 Re-reference (average)
EEG = pop_reref( EEG, []);
         
% Line#15 Saving dataset before extracting events
EEG = pop_saveset( EEG, 'filename','xxx.set','filepath','C:\\folder\\location\\Subject_01\\');
