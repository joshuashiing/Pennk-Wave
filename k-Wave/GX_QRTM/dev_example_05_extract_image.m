clear;clc

migdir = 'Data_example_05/';
mig_list = {'mig_cd_rf.mat', 'mig_cd_ac.mat', 'mig_cd.mat', 'mig_cd_sq.mat'};
% mig_list = {'mig_cd_rf.mat', 'mig_cd_ac.mat'};

% reference
mig_file = [migdir 'mig_cd_rf.mat'];
d = load(mig_file);
mig1_rf = del2(sum(d.mig1, 3));
mig2_rf = del2(sum(d.mig2, 3));

% acoustic
mig_file = [migdir 'mig_cd_ac.mat'];
d = load(mig_file);
mig1_ac = del2(sum(d.mig1, 3));
mig2_ac = del2(sum(d.mig2, 3));

% visco-acoustic
mig_file = [migdir 'mig_cd.mat'];
d = load(mig_file);
mig1 = del2(sum(d.mig1, 3));
mig2 = del2(sum(d.mig2, 3));

% visco-acoustic smooth Q
mig_file = [migdir 'mig_cd_sq.mat'];
d = load(mig_file);
mig1_sq = del2(sum(d.mig1, 3));
mig2_sq = del2(sum(d.mig2, 3));

save([migdir 'rtm_images.mat'], 'mig1_rf', 'mig2_rf', 'mig1_ac',...
    'mig2_ac', 'mig1', 'mig2', 'mig1_sq', 'mig2_sq');