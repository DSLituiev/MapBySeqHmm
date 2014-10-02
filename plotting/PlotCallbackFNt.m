function output_txt = PlotCallbackFNt(~, event_obj, lser, AllReads,  xname, names, flaglog, selection)
% Display the position of the data cursor
% obj          Currently not used (empty)
% event_obj    Handle to event object
% output_txt   Data cursor text string (string or cell array of strings).
% lser         Line series handle vector
% flaglog      boolean: whether to take the log10 of the y-data

if isempty(selection)
    selection = true(numel(AllReads.chromosome), 1);
end

dataSubs = get(event_obj,'DataIndex');

% pos = get(event_obj,'Position');
targ = event_obj.Target;

[chr, figN] = find(lser==targ, 1, 'first');

if isempty(chr)
    output_txt = '';
    return
end

output_txt = {sprintf('%s;\t', names{figN})};

ChrReads.(xname) = AllReads.(xname)(AllReads.chromosome == chr & selection);

ChrReads.(names{figN}) = AllReads.(names{figN})(AllReads.chromosome == chr & selection);
if isprop(AllReads, 'r') && isprop(AllReads, 'f')
    ChrReads.f = AllReads.f(AllReads.chromosome == chr & selection);
    ChrReads.r = AllReads.r(AllReads.chromosome == chr & selection);
end


if isprop(AllReads, 'geneID')
    ChrReads.geneID = AllReads.geneID(AllReads.chromosome == chr & selection);
    ChrReads.geneSO = AllReads.geneSO(AllReads.chromosome == chr & selection);
end
indAllR = find( AllReads.chromosome == chr & selection & AllReads.(xname) == ChrReads.(xname)(dataSubs));

% fprintf('x:\t%8.0f, P:\t%2.2f %%\n', [round(pos(1).*1e6), pos(2).*100])
if flaglog
    
    fprintf('chr:\t%u\t',  chr)
    fprintf([xname, ':\t%8.0f\t'], ChrReads.(xname)(dataSubs))
    fprintf([names{figN}, ':\t%5.3f\t'], ChrReads.(names{figN})(dataSubs) )
    fprintf(['index', ':\t%5.0u\t'], indAllR )
    output_txt = { output_txt{:}, strcat(sprintf('X\t: %u ,  \tlog10(Y)\t: %4.3f',  ChrReads.(xname)(dataSubs), (ChrReads.(names{figN})(dataSubs)))) };     
    
else
    fprintf('chr:\t%u\t',  chr)
    fprintf([xname, ':\t%8.0f\t'], ChrReads.(xname)(dataSubs))
    fprintf([names{figN}, ':\t%5.3f\t'], ChrReads.(names{figN})(dataSubs) )
    fprintf(['index', ':\t%5.0u\t'], indAllR )
    output_txt = {output_txt{:}, sprintf('X\t: %u  ,  \tY\t: %4.3f',...
        ChrReads.(xname)(dataSubs), ChrReads.(names{figN})(dataSubs) )};  
    
end

if isfield(ChrReads, {'r', 'f'})
    fprintf(['SNP-ratio (f)', ':\t%4.3f\t'], ChrReads.f(dataSubs) )
    output_txt = { output_txt{:},...
        sprintf('f\t: %4.3f ,  \tr\t: %u', ChrReads.f(dataSubs), ChrReads.r(dataSubs))};
end

if isfield(ChrReads, 'geneID')
    fprintf(['geneID', ':\t%s\t'], ChrReads.geneID{dataSubs} )
    fprintf(['geneSO', ':\t%u\t'], ChrReads.geneSO(dataSubs) )
    output_txt = { output_txt{:},...
        sprintf('EffectSO: %u', ChrReads.geneSO(dataSubs)),...
        sprintf(' geneID : %s', ChrReads.geneID{dataSubs})};
end

fprintf('\n')
% % If there is a Z-coordinate in the position, display it as well
% if length(pos) > 2
%     output_txt{end+1} = ['Z: ',num2str(pos(3),4)];
% end
