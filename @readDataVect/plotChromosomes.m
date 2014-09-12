function  plotChromosomes(obj, fieldName, varargin)
            p = inputParser;
            addRequired(p, 'obj',@isobject);
            addRequired(p, 'fields', @ischar);
            %
            addParamValue(p,     'exp10',             false, @isscalar);
            % addParamValue(p,     'color',                   'b');
            addParamValue(p,     'xname',                  'x',  @ischar);
            addParamValue(p,     'xlim',   0,  @isnumeric);
            addParamValue(p,   'FigName',                   '',  @ischar);
            addParamValue(p,     'figure',  'old',  @ischar);
            addParamValue(p,    'OldFig',   zeros(obj.chrNumber,1),  @isnumeric);
            addParamValue(p,    'yscale',               'lin' ,  @ischar);
            addParamValue(p,    'xscale',               'lin' ,  @ischar);
            addParamValue(p,      'ylim',                   []);
            addParamValue(p, 'linestyle',                  '-',  @ischar);
            addParamValue(p,   'plotfun',                @plot,  @(x)isa(x,'function_handle'));
            addParamValue(p,    'select',                   '',  @ischar );
            addParamValue(p,   'linewidth',                1 , @isscalar);
            addParamValue(p,    'norm',                false , @isscalar);
            addParamValue(p,    'yThr',                NaN , @isscalar);
            
            parse(p, obj, fieldName, varargin{:});
            
            if isempty(obj.(fieldName))
               warning('plotChromosomes:emptyInputField', 'the input field is empty!')
               return
            end
            
            if isempty(p.Results.FigName)
                if iscellstr(p.Results.fields)
                    if strcmp(p.Results.xname,'x')
                        FigureName = strjoin( p.Results.fields, ', ' );
                    else
                        FigureName = [ strjoin( p.Results.fields, ', ' ), ' vs ', p.Results.xname];
                    end
                elseif  ischar(p.Results.fields)
                    if strcmp(p.Results.xname,'x')
                        FigureName = p.Results.fields;
                    else
                        FigureName = [ p.Results.fields, ' vs ', p.Results.xname];
                    end
                end
            else
                FigureName = p.Results.FigName;
            end
            flagFigOld = strcmpi(p.Results.figure, 'old') && ~isempty(obj.prevFig) && ...
                ishghandle(obj.prevFig(end)) && all(ishghandle(obj.prevAxes(:,end)));% && isvalid(obj.prevFig) ;
            
            if  flagFigOld          
                f = figure(obj.prevFig(end));
                set(f, 'name', FigureName);
                %  Params.Colors = ['b','g'];
            else
                f = figure('name', FigureName);
                % Params.Colors = ['y','r'];
            end

            maxNt = zeros(obj.chrNumber, 1);
            maxXMb = zeros(obj.chrNumber, 1);
            spl = zeros(obj.chrNumber, 1);
            
            nPlLines = numel(obj.prevYNames);
            
            lserClmn = nPlLines + 1; % size(obj.prevLine, 2)+1;
            obj.prevYNames{end+1} = fieldName;
            
            %% plot
            maxY = -Inf;
            minY = Inf;
            for chr = 1:obj.chrNumber
                inds = (obj.chromosome == chr); 
                
                if flagFigOld
                    spl(chr) =   obj.prevAxes(chr, end);
                    subplot(spl(chr));
                else
                    spl(chr) = subplot(double(obj.chrNumber), 1, double(chr));
                end
                
                if isscalar(obj.(fieldName))
                    xValues = obj.(p.Results.xname)(inds);
                    xValues = xValues([1,end]);
                    yValues = [1, 1] * obj.(fieldName);
                elseif isempty(p.Results.select)
                    xValues = obj.(p.Results.xname)(inds);
                    yValues = obj.(fieldName)(inds);
                else
                    xValues = obj.(p.Results.xname)(p.Results.select);
                    yValues = obj.(fieldName)(p.Results.select);
                end
               
                
                if p.Results.norm
                    yValues = yValues + obj.cNormConst(chr);
                end
                
                maxNt(chr) = max( xValues );
                maxXMb(chr) = single(maxNt(chr))*1e-6;
                
                if ( numel(xValues) == numel(yValues) )
                    obj.prevLine(chr, lserClmn) = doplottingcurve( xValues, yValues, inds, p );
                elseif isscalar(yValues)
                    obj.prevLine(chr, lserClmn) = doplottingline( maxNt(chr), yValues, inds, p );
                else
                    warning('plotAllChr:dimmismatch','dimension mismatch between .x and the .(%s)! Nothing to plot', fieldName)
                end
                
                hold all
                
                if ~isnan(p.Results.yThr)
                    plot( xValues([1,end]), p.Results.yThr.*[1,1], 'k--')
                end
                    
                maxY = nanmax([maxY; p.Results.yThr; yValues(:)]);
                minY = nanmin([minY; p.Results.yThr; yValues(:)]);
                
                %                 spl(chr) = gca;
                if strcmpi(p.Results.xscale , 'lin')
                    set(spl(chr), 'xlim', [ 0 , maxNt(chr) ]);
                else
                    %  set(spl(chr),'xlim',  10.^[0 , log10( double(ChrReads(chr).maxNt) ) ]);
                    set(spl(chr), 'xlim', 10.^[0, round(log10( double(maxNt(chr)) ))] )
                end
            end % loop through the chromosomes
                        
            set(obj.prevLine(:, lserClmn), 'linewidth', p.Results.linewidth);
            
            set(spl,'yscale',p.Results.yscale);
            set(spl,'xscale',p.Results.xscale);
            
            minY = floor(minY);
            maxY = ceil(maxY);
            
            if ~isempty(p.Results.ylim) && isnumeric(p.Results.ylim)
                obj.prevYLims = p.Results.ylim;
            elseif  ~flagFigOld % isempty(p.Results.ylim)
                obj.prevYLims = [minY, maxY];
            else
                obj.prevYLims(1) = nanmin(obj.prevYLims(1), minY);
                obj.prevYLims(2) = nanmax(obj.prevYLims(2), maxY);
            end

            set(spl, 'ylim', obj.prevYLims);
            %    set(spl,'YTickMode','auto')
            if strcmpi(p.Results.yscale, 'log')
                if p.Results.exp10
                    set(spl, 'ytick',...
                        10.^( (obj.prevYLims(1)) : floor(diff((obj.prevYLims))/2) : (obj.prevYLims(2)) ) );
                else
                    set(spl, 'ytick',...
                        10.^( log10(obj.prevYLims(1)) : floor(diff(log10(obj.prevYLims))/2) : log10(obj.prevYLims(2)) ) );
                end
            end
            
            
            %% set the x-length in correspondence if ( p.Results.xname == 'x' )
            %        [left bottom width height]
            if ~isscalar(spl)
                splpos = cell2mat(get(spl,'position'));
            else
                splpos = get(spl,'position');
            end
            if ~(size(p.Results.xlim,2) == 2)
                if strcmpi(p.Results.xname,'x')
                    splpos(:,3) = max(splpos(:,3)).* ([maxXMb(:)])'./max([maxXMb(:)]);
                    % elseif strcmpi(p.Results.xscale, 'lin')
                    %     splpos(:,3) = splpos(1,3).* ([ChrReads(:).maxXMb])'./max([ChrReads(:).maxXMb]);
                elseif strcmpi(p.Results.xscale, 'log')
                    % xlims = cellfun(@diff, get(spl, 'xlim'))
                    splpos(:,3) = max(splpos(:,3)).*...
                        ( log10(double(maxNt'))./ log10(double(max(maxNt))) );
                end
            end
            
            if strcmpi(p.Results.xname,'x')
                xticks = 1e6*(0:1:ceil(max(maxXMb)))';
                xticklabel =  strtrim(cellstr(num2str(1e-6 * xticks,3)));
                ind5 = logical(mod(xticks,5e6));
                xticklabel(ind5) = {''};
                
                for chr = 1:obj.chrNumber
                    % set(spl(chr), 'position', splpos(chr,:));
                    % get(spl(chr),'position')
                    set(spl(chr),'xtick',xticks(1:ceil( maxXMb(chr))),...
                        'XTickLabel', xticklabel(1:ceil( maxXMb(chr))));
                    set(spl(chr),'box','off')
                end
            end
            
            if (size(p.Results.xlim,2) == 2)
                set(spl, 'xlim', p.Results.xlim);
                for chr = 1:obj.chrNumber
                    splpos(chr,3) = max(splpos(:,3));
                    set(spl(chr), 'position', splpos(chr,:) );
                end
            else
                for chr = 1:obj.chrNumber
                    set(spl(chr), 'position', splpos(chr,:));
                end
            end
            %== define the callback function:
            dcm_obj = datacursormode(f);
            set(dcm_obj ,'UpdateFcn',{@PlotCallbackFNt, obj.prevLine, obj, p.Results.xname, obj.prevYNames, p.Results.exp10});
            
            
            %% subfunctions
            function lser = doplottingcurve(xValues , yValues , inds, p)
                if nargin(p.Results.plotfun) < 3
                    if    ~p.Results.exp10
                        lser = p.Results.plotfun(xValues, yValues);
                    elseif p.Results.exp10
                        lser = p.Results.plotfun( xValues, 10.^yValues);
                    end
                else
                    if    ~p.Results.exp10
                        lser = p.Results.plotfun(xValues, yValues, inds);
                    elseif p.Results.exp10
                        lser = p.Results.plotfun( xValues, 10.^yValues, inds);
                    end
                end
            end
            
            function lserOut = doplottingline(maxNt, yValues,  p)
                if    ~p.Results.exp10
                    lserOut =  p.Results.plotfun([0 1] .*double(maxNt), yValues .*[1 1]);
                elseif p.Results.exp10
                    lserOut =  p.Results.plotfun([0 1] .*double(maxNt), 10.^double(yValues) .*[1 1]);
                end
            end
        
        if ~flagFigOld
                  obj.prevFig(end+1) = f;
        end
                  obj.prevAxes(:, nPlLines+1) = spl;
        end