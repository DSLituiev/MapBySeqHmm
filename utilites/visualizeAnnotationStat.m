function visualizeAnnotationStat(annotation)

F = fields(annotation);
if isempty(F)
   warning('visualizeAnnotationStat:emptyAnnotation', 'the annotation is empty!') 
else
    
for ii = numel(F):-1:1
%     S.(F{ii}) = sum( cellfun( @(x)~isempty(x), annotation.(F{ii})) );
    numF(ii) =  sum( cellfun( @(x)~isempty(x), annotation.(F{ii})(:,1) ) );
    nameF{ii} = F{ii}(4:end);
    var0 = regexp(nameF{ii}, '_variant', 'split');
    nameF{ii} = strcat(var0{:}, sprintf(': %u', numF(ii)));
end

[numFs, ind] = sort(numF, 'descend');
nameFs = nameF(ind);
figure
barh( numFs(logical(numFs)) )
set(gca, 'yticklabel', nameFs(logical(numFs)))
set(gca, 'position', [0.27, 0.1, .67, .8])
set(gca, 'xscale', 'log')
% set(hx, 'interpreter', 'none')

end