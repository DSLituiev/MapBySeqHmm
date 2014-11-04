%                               431-400 Year Long Project 
%                               LA1 - Medical Image Processing 2003
% Supervisor     :  Dr Lachlan Andrew
% Group Members  :  Alister Fong    78629   a.fong1@ugrad.unimelb.edu.au
%                   Lee Siew Teng   102519  s.lee1@ugrad.unimelb.edu.au
%                   Loh Jien Mei    103650  j.loh1@ugrad.unimelb.edu.au
%
% File and function name : finishing_beep
% Version                : 1.0
% Date of completion     : 18 July 2003   
% Written by    :   Alister Fong    78629   a.fong1@ugrad.unimelb.edu.au
%
% Input   :
%         (ALL OPTIONAL, EITHER ENTER THEM IN ORDER OR DON'T ENTER AT ALL)
%         ------------------------------
%         PARAMETER VARIABLE SEQUENCE 1
%         ------------------------------
%         sequence    -   Sequence number to be played in turn 
% 						sequence choices
%                             1 - train
%                             2 - gong
%                             3 - splat
%                             4 - handel
%                             5 - laughter
%                             6 - chirp
%                             7 - beeps
%         pause_time  -   Pause time before the next sequence is played
%         beep_times  -   Number of beeps if sequence 7 (beeps) is selected
%         
%         All these matrices if entered, must be of the same length
% 			
%                 OR
%         ------------------------------
%         PARAMETER VARIABLE SEQUENCE 2
%         ------------------------------
% 			'playset',
% 			playset_number -
% 								1 - Mixed combination of sequences
% 								2 - Beeping sequence
% 								3 - Continuous gong sounds
% 								4 - Continuous train sounds
% 								5 - Laughter sequence
% 								6 - splat and Laughter sequence
%
% 
%         DEFAULT:    Random selection of preset sequences
%
% Description:
%     Just some light hearted sounds intended to signify the end of a program
%     
% Usage >> finishing_beep or
%          finishing_beep(sequence,pause_time,beep_times) or
% 	 	   finishing_beep('playset',playset_number)

function finishing_beep(varargin)
% Default leeway for each of the standard sounds to allow for loading time
leeway = 0.3;
% Default setting
playset    = -1;
sequence   = [];
pause_time = [];
beep_times = [];

if ~isempty(varargin)
    if length(varargin) == 3
        playset    = 0;
        sequence   = varargin{1};
        pause_time = varargin{2};
        beep_times = varargin{3};
    elseif length(varargin) == 2
        if strcmp(varargin{1},'playset') == 1
            playset = varargin{2};
        else
            error('Invalid parameters');
        end
    else
        error('Invalid number of variables');
    end
end
if  playset ~= 0
    if playset < 0
        playset = randperm(5);
    end
	switch playset(1)
	case 1
        % Mixed combination of sequences
        sequence   = [1,2,3,4,5];
		pause_time = [0,0,0,0,0];
		beep_times = [0,0,0,0,0];
	case 2
        % Beeping sequence
		pause_time = [0.1:0.1:1];
		beep_times = [10:10:100];
        sequence   = pause_time;
        sequence(:)= 7;
    case 3
        % Continuous gong sounds
        beep_times = [1:1:10];
        pause_time = beep_times;
        pause_time(:)=0;
        sequence   = pause_time;
        sequence(:)= 2;
    case 4
        % Continuous train sounds
        beep_times = [1:1:10];
        pause_time = beep_times;
        pause_time(:)=0;
        sequence   = pause_time;
        sequence(:)= 1;
    case 5
        % Laughter sequence
        beep_times = [1:1:10];
        pause_time = beep_times;
        pause_time(:)=0;
        sequence   = pause_time;
        sequence(:)= 5;
    case 6
        % splat and Laughter sequence
        beep_times = [0,0];
        sequence   = [3,5];
        pause_time = [0,0];
	otherwise
        error('Invalid playset selection');
	end
end

% -------------------------------------------
% Check the length of the matrices
% -------------------------------------------
if (length(sequence) ~= length(pause_time)) || (length(sequence) ~= length(beep_times))
    error('sequence length not equal to pause_time length');
end

% -------------------------------------------
% Format the matrices properly
% -------------------------------------------
sequence = [-1,sequence];
sequence_length = length(sequence);
pause_time = [-1,pause_time];
beep_times = [-1,beep_times];

% -------------------------------------------
% Display message
% -------------------------------------------
% disp('Press Ctrl-C to stop otherwise wait for sound sequence to end');

% -------------------------------------------
% Play sequences and delays
% -------------------------------------------
for c = 2:1:sequence_length
    try
		switch sequence(c)
		case 1  %train
            if sequence(c-1) ~= sequence(c)
                load train;
            end
            pause_time(c) = pause_time(c) + 1.5 + leeway;
		case 2 %gong
            if sequence(c-1) ~= sequence(c)
                load gong;
            end
            pause_time(c) = pause_time(c) + 5.1 + leeway;
		case 3 %splat
            if sequence(c-1) ~= sequence(c)
                load splat
            end
            pause_time(c) = pause_time(c) + 1.2 + leeway;
		case 4 %handel
            if sequence(c-1) ~= sequence(c)
				load handel;
            end
            pause_time(c) = pause_time(c) + 8.9 + leeway;
		case 5 %laughter
            if sequence(c-1) ~= sequence(c)
        		load laughter
            end
            pause_time(c) = pause_time(c) + 6.4 + leeway;
		case 6 %chirps
			% There is something wrong with this one sometimes, won't stop!!
            if sequence(c-1) ~= sequence(c)
				load chirp;
            end
            pause_time(c) = pause_time(c) + 1.6 + leeway;
        case 7 %beeps
            for i=1:1:beep_times(c)
                beep;
            end
        otherwise
            error('Invalid sequence choice');
		end
		
        %Play the sounds except for beeps    
        if sequence(c) ~= 7
            wavplay(y,Fs,'async');
        end

        %pause time before the next sequence
        pause(pause_time(c));    
    catch
		% If the user presses 'Ctrl-C' to terminate then we do not want an error message
        % But want to exit from this program
        return;
    end
end