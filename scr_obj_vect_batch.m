close all; clear all; clc;
dbclear if warning
tic


DATA_PATH = '/media/Processing/seq/data';
% DATA_PATH = '/media/Processing/seq/olddata';
% DATA_PATH = './data';
USERFNCT_PATH = '/media/Processing/MATLABuserfunctions';
% addpath(fullfile(USERFNCT_PATH, 'MATLABuserfunctions/binomial') );
addpath(USERFNCT_PATH);
addpath(fullfile(USERFNCT_PATH, 'mtimesx'));
addpath(fullfile(USERFNCT_PATH, 'MinMaxSelection'));
addpath(fullfile(USERFNCT_PATH, 'newtonraphson'));
addpath(fullfile(USERFNCT_PATH, 'fastmedfilt1d'));
savepath;
addpath('./utilites');
addpath('./plotting');
addpath('./emission');

%=  provide the reference ID, which is the name of the csv file without '.csv'
%= together with the backround genotype ID
%= known positions of the causative SNP can be also provided here for
%= further visualization

%= number of plants:

% dataID = 'CH/73-bwa-rmdup-clipOverlap-q20-nm2-ems-annotation'; N = 100;
% dataID = 'CH/91-bwa-rmdup-clipOverlap-q20-nm2-ems-annotation'; N = 100;
% dataID = 'CH/91-ngm-rmdup-clipOverlap-q20-ems-annotation'; N = 100;
%  dataID = 'CH/39-bwa-rmdup-clipOverlap-q20-nm2-ems-annotation'; N = 100;
% dataID = 'CH/39-bwa-rmdup-clipOverlap-q20-ems-annotation'; N = 100;

% dataID = '91-ngm-rmdup-clipOverlap-q20-ems-annotation'; N = 100;
% dataID = 'ABD/20140605.X-ABD192-bwa-rmdup-clipOverlap-q20-nm2-ems-annotation';  N = 100;

% dataList = {'ABD/20120427-ABD159-ngm-rmdup-clipOverlap-q20-nm2-ems-annotation',   100,  2,  17521246 ;...
%             'ABD/20120427-ABD159-bwa-rmdup-clipOverlap-q20-nm2-ems-annotation',   100,  2,  17521246 ;...
%             'ABD/20130426.B-ABD173-ngm-rmdup-clipOverlap-q20-nm2-ems-annotation',  50,  3, 1619248; ...
%             'ABD/20130426.B-ABD173-bwa-rmdup-clipOverlap-q20-nm2-ems-annotation',  50 , 3, 1619248; ...
%             'HL10/HL10-rmdup-clipOverlap-nm2-q20-temp-ems-annotation', 50, 3 , 16473265; ...
%             'HL10/HL10-rmdup-clipOverlap-nm2-q20-ems-annotation',  50, 3 , 16473265; ...
%             'HL10/HL10-rmdup-clipOverlap-q20-freebayes-ems-annotation',  50, 3 , 16473265; ...
%             'HL7/p889_20110125_HL7_Paired-rmdup-clipOverlap-nm2-q20-temp-ems-annotation', 50, 1,  5672441; 
%             'HL7/p889_20110125_HL7_Paired-rmdup-clipOverlap-nm2-q20-ems-annotation', 50, 1,  5672441; 
%             'HL7/p889_20110125_HL7_Paired-rmdup-clipOverlap-q20-freebayes-ems-annotation', 50, 1, 5672441};
            
           
% dataID = 'ABD/20120427-ABD159-ngm-rmdup-clipOverlap-q20-nm2-ems-annotation';   N = 100; chr0 = 2; x0 = 17521246;
% dataID = 'ABD/20120427-ABD159-bwa-rmdup-clipOverlap-q20-nm2-ems-annotation'; chr0 = 2; x0 = 17521246;  N = 100;

%  dataID = 'ABD/20130426.B-ABD173-ngm-rmdup-clipOverlap-q20-nm2-ems-annotation'; chr0 = 3 ; x0 = 1619248;  N = 50;  % bkgrID = 'ABD241-rmdup-clipOverlap-q20-freebayes';
%  dataID = 'ABD/20130426.B-ABD173-bwa-rmdup-clipOverlap-q20-nm2-ems-annotation'; chr0 = 3 ; x0 = 1619248;  N = 50;  % bkgrID = 'ABD241-rmdup-clipOverlap-q20-freebayes';
 
% dataID = 'MO8_mutant_pool-rmdup-clipOverlap-q20-freebayes-ems-annotation-repfilt';

% dataID = 'HL4/p889_20110501_HL4_pairing-rmdup-clipOverlap-nm2-q20-ems-annotation';  chr0 =  3 ;  x0 =  16473265;N = 50;
% dataID = 'HL10/HL10-rmdup-clipOverlap-nm2-q20-temp-ems-annotation';  chr0 =  3 ;  x0 =  16473265;N = 50;
% dataID = 'HL10/HL10-rmdup-clipOverlap-nm2-q20-ems-annotation';  chr0 =  3 ;  x0 =  16473265;N = 50;
% dataID = 'HL10/HL10-rmdup-clipOverlap-q20-freebayes-ems-annotation';  chr0 =  3 ;  x0 =  16473265;N = 50;

% dataID = 'HL7/p889_20110125_HL7_Paired-rmdup-clipOverlap-nm2-q20-temp-ems-annotation'; x0 = 5672441; chr0 = 1; N = 50;
% dataID = 'HL7/p889_20110125_HL7_Paired-rmdup-clipOverlap-nm2-q20-ems-annotation'; x0 = 5672441; chr0 = 1; N = 50;
% dataID = 'HL7/p889_20110125_HL7_Paired-rmdup-clipOverlap-q20-freebayes-ems-annotation'; x0 = 5672441; chr0 = 1; N = 50;

dataList = {'CH/73-bwa-rmdup-clipOverlap-q20-nm2-ems-annotation', 100, [], []};
    
ALPHA = 4;
% AR.Alpha = 1./(0:0.01:1);

for ii = 1:size(dataList,1)
    run_obj_hmm(dataList(ii,:), ALPHA )
end

%%
