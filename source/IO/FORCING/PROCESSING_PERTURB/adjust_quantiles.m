%========================================================================
% CryoGrid FORCING post-processing class set_min_windspeed
%
%
% Authors:
% S. Westermann, December 2022
%
%========================================================================

classdef adjust_quantiles < matlab.mixin.Copyable 
    
    properties
        CONST
        PARA
        STATVAR
        TEMP
    end
    
    methods
        function proc = provide_PARA(proc)
            
            proc.PARA.variable = [];
            proc.PARA.number_of_bins = [];
            proc.PARA.method = []; %histogram or quantiles
            proc.PARA.bin_fraction = [];
            proc.PARA.bin_adjustment = [];                   
        end
        
        
        function proc = provide_CONST(proc)

        end
        
        
        function proc = provide_STATVAR(proc)
            
        end
        
        
        function proc = finalize_init(proc, tile)

        end
        
        
        function forcing = process(proc, forcing, tile)

            if strcmp(proc.PARA.method, 'histogram')
                [N,edges,bin] = histcounts(forcing.DATA.(proc.PARA.variable), proc.PARA.number_of_bins);
            elseif strcmp(proc.PARA.method, 'quantiles')
                edges = [-Inf quantile(forcing.DATA.(proc.PARA.variable), proc.PARA.number_of_bins-1) Inf]';
                bin = forcing.DATA.(proc.PARA.variable) .* NaN;
                for j=1:proc.PARA.number_of_bins
                    bin(find(forcing.DATA.(proc.PARA.variable)>= edges(j) & forcing.DATA.(proc.PARA.variable) < edges(j+1))) = j;
                end
            end
                
            proc.TEMP.bin = bin;
            
            factor = linspace(0,1,proc.PARA.number_of_bins);
            factor2 = factor;
            for i=2:size(proc.PARA.bin_adjustment,1)
                
                factor2(factor>= proc.PARA.bin_fraction(i-1) & factor<=proc.PARA.bin_fraction(i)) = proc.PARA.bin_adjustment(i-1) + ...
                    (proc.PARA.bin_adjustment(i) - proc.PARA.bin_adjustment(i-1)) .* (factor(factor>= proc.PARA.bin_fraction(i-1) & factor<=proc.PARA.bin_fraction(i)) - proc.PARA.bin_fraction(i-1));
            end
            
            for i=1:proc.PARA.number_of_bins
                forcing.DATA.wind(bin==i) = forcing.DATA.wind(bin==i) .* factor2(i);
            end


        end
        
        %----------------------------------------------
        %perturb
        function proc = preprocess(proc, forcing, tile)

            if strcmp(proc.PARA.method, 'histogram')
                [N,edges,bin] = histcounts(forcing.DATA.(proc.PARA.variable), proc.PARA.number_of_bins);
            elseif strcmp(proc.PARA.method, 'quantiles')
                edges = [-Inf quantile(forcing.DATA.(proc.PARA.variable), proc.PARA.number_of_bins-1) Inf]';
                bin = forcing.DATA.(proc.PARA.variable) .* NaN;
                for j=1:proc.PARA.number_of_bins
                    bin(find(forcing.DATA.(proc.PARA.variable)>= edges(j) & forcing.DATA.(proc.PARA.variable) < edges(j+1))) = j;
                end
            end
                
            proc.TEMP.bin = bin;
            
            factor = linspace(0,1,proc.PARA.number_of_bins);
            factor2 = factor;
            for i=2:size(proc.PARA.bin_adjustment,1)
                
                factor2(factor>= proc.PARA.bin_fraction(i-1) & factor<=proc.PARA.bin_fraction(i)) = proc.PARA.bin_adjustment(i-1) + ...
                    (proc.PARA.bin_adjustment(i) - proc.PARA.bin_adjustment(i-1)) .* (factor(factor>= proc.PARA.bin_fraction(i-1) & factor<=proc.PARA.bin_fraction(i)) - proc.PARA.bin_fraction(i-1));
            end
            
            proc.STATVAR.bins = edges;
            proc.STATVAR.adjustment_factor = factor2';

        end
        
        function forcing = perturb_forcing(proc, forcing, tile)
            
            forcing.TEMP.wind = forcing.TEMP.wind .* proc.STATVAR.adjustment_factor(find(forcing.TEMP.wind - proc.STATVAR.bins <= 0, 1, 'first')-1,1);

        end
        
        
                %-------------param file generation-----
%         function proc = param_file_info(post_proc)
%             post_proc = provide_PARA(post_proc);
% 
%             post_proc.PARA.STATVAR = [];
%             post_proc.PARA.class_category = 'FORCING POST_PROCESSING';
%             post_proc.PARA.options = [];
%             
%             post_proc.PARA.eliminate_fraction = [];
%             post_proc.PARA.survive_fraction = [];
%                         
%             post_proc.PARA.default_value.window_size = {7};
%             post_proc.PARA.comment.window_size = {'window size in days within which precipitation is reallocated'};
%             
%             post_proc.PARA.default_value.eliminate_fraction = {0.5};
%             post_proc.PARA.comment.eliminate_fraction = {'fraction of smallest precipitation events (= timestamps with precipitation) that is reallocated to larger events'};
%             
%             post_proc.PARA.default_value.survive_fraction = {0.5};  
%             post_proc.PARA.comment.survive_fraction = {'fraction of largest precipitation events (= timestamps with precipitation) that the small events are reallocated to'};
%             
%         end
        
    end
    
end

