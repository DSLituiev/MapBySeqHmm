classdef soDictionary
    properties
        SO
        descr
    end
    methods
        function obj = soDictionary(varargin)
            if nargin>0
                soFilePath = varargin{1};
            else
                soFilePath = './reference/SO_terms.csv';
            end
            fid = fopen(soFilePath);
            assert(fid>0, sprintf('the SO file was not found. the specified location was: \n%s\n', soFilePath) )
            terms = textscan(fid, '%s %s %s %s %s %s', 'delimiter', ';', 'HeaderLines', 0);
            fclose(fid);
            a = char(terms{3}{:});
            obj.SO = str2num(a(:,4:end));
            obj.descr = terms{1};
            clear a terms
        end
        function currEffect = lookup(obj, inputSO)
             currEffect = obj.descr{obj.SO ==inputSO};
        end
    end
end