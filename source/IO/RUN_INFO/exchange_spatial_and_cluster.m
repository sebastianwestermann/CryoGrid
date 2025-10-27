%========================================================================
% CryoGrid exchange class exchange_spatial_and_cluster

% S. Westermann, April 2025
%========================================================================

classdef exchange_spatial_and_cluster < matlab.mixin.Copyable
    
    properties
        PARA
        CONST
        STATVAR
    end
    
    
    methods
        
        
        function spatial_exchange = provide_PARA(spatial_exchange)
            
            spatial_exchange.PARA.read_variables = [];
            spatial_exchange.PARA.write_variables = [];
            
        end
        
        function spatial_exchange = provide_CONST(spatial_exchange)
            
        end
        
        function spatial_exchange = provide_STATVAR(spatial_exchange)
            
        end
        
        
        function spatial_exchange = finalize_init(spatial_exchange, run_info)

        end
        

        function run_info = write_spatial_to_master(spatial_exchange, current_run_info, run_info)
            if ~isempty(spatial_exchange.PARA.write_variables) && iscell(spatial_exchange.PARA.write_variables)
                for i=1:size(spatial_exchange.PARA.write_variables,1)
                     run_info.STORE.(spatial_exchange.PARA.write_variables{i,1}) = current_run_info.(spatial_exchange.PARA.write_variables{i,1});
                end
            end
        end
        

        function current_run_info = read_spatial_from_master(spatial_exchange, current_run_info, run_info)
            if ~isempty(spatial_exchange.PARA.read_variables) && iscell(spatial_exchange.PARA.read_variables)
                for i=1:size(spatial_exchange.PARA.read_variables,1)
                    current_run_info.(spatial_exchange.PARA.read_variables{i,1}) = run_info.STORE.(spatial_exchange.PARA.read_variables{i,1});
                end
            end
        end
        
    end
end



