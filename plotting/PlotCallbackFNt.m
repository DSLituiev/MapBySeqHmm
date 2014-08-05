function output_txt = PlotCallbackFNt(~, event_obj, lser, AllReads,  xname, names, flaglog)
% Display the position of the data cursor
% obj          Currently not used (empty)
% event_obj    Handle to event object
% output_txt   Data cursor text string (string or cell array of strings).
% lser         Line series handle vector
% flaglog      boolean: whether to take the log10 of the y-data

dataIndex = get(event_obj,'DataIndex');
% pos = get(event_obj,'Position');
targ = event_obj.Target;

[chr, figN] = find(lser==targ, 1, 'first');

if isempty(chr)
    output_txt = '';
    return
end

output_txt = {sprintf('%s;\t', names{figN})};

ChrReads.(xname) = AllReads.(xname)(AllReads.chromosome == chr);
ChrReads.(names{figN}) = AllReads.(names{figN})(AllReads.chromosome == chr);
if isprop(AllReads, 'r') && isprop(AllReads, 'f')
    ChrReads.f = AllReads.f(AllReads.chromosome == chr);
    ChrReads.r = AllReads.r(AllReads.chromosome == chr);
end


if isprop(AllReads, 'geneID')
    ChrReads.geneID = AllReads.geneID(AllReads.chromosome == chr);
    ChrReads.geneSO = AllReads.geneSO(AllReads.chromosome == chr);
end
indAllR = find( AllReads.chromosome == chr & AllReads.(xname) == ChrReads.(xname)(dataIndex));

% fprintf('x:\t%8.0f, P:\t%2.2f %%\n', [round(pos(1).*1e6), pos(2).*100])
if flaglog
    
    fprintf('chr:\t%u\t',  chr)
    fprintf([xname, ':\t%8.0f\t'], ChrReads.(xname)(dataIndex))
    fprintf([names{figN}, ':\t%5.3f\t'], ChrReads.(names{figN})(dataIndex) )
    fprintf(['index', ':\t%5.0u\t'], indAllR )
    output_txt = { output_txt{:}, strcat(sprintf('X\t: %u ,  \tlog10(Y)\t: %4.3f',  ChrReads.(xname)(dataIndex), (ChrReads.(names{figN})(dataIndex)))) };     
    
else
    fprintf('chr:\t%u\t',  chr)
    fprintf([xname, ':\t%8.0f\t'], ChrReads.(xname)(dataIndex))
    fprintf([names{figN}, ':\t%5.3f\t'], ChrReads.(names{figN})(dataIndex) )
    fprintf(['index', ':\t%5.0u\t'], indAllR )
    output_txt = {output_txt{:}, sprintf('X\t: %u  ,  \tY\t: %4.3f',...
        ChrReads.(xname)(dataIndex), ChrReads.(names{figN})(dataIndex) )};  
    
end

if isfield(ChrReads, {'r', 'f'})
    fprintf(['SNP-ratio (f)', ':\t%4.3f\t'], ChrReads.f(dataIndex) )
    output_txt = { output_txt{:},...
        sprintf('f\t: %4.3f ,  \tr\t: %u', ChrReads.f(dataIndex), ChrReads.r(dataIndex))};
end

if isfield(ChrReads, 'geneID')
    fprintf(['geneID', ':\t%s\t'], ChrReads.geneID{dataIndex} )
    fprintf(['geneSO', ':\t%u\t'], ChrReads.geneSO(dataIndex) )
    output_txt = { output_txt{:},...
        sprintf('EffectSO: %u', ChrReads.geneSO(dataIndex)),...
        sprintf(' geneID : %s', ChrReads.geneID{dataIndex})};
end

fprintf('\n')
% % If there is a Z-coordinate in the position, display it as well
% if length(pos) > 2
%     output_txt{end+1} = ['Z: ',num2str(pos(3),4)];
% end
