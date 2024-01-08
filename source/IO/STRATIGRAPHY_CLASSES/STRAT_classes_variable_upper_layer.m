%========================================================================
% CryoGrid STRATIGRAPHY_CLASSES class STRAT_classes 
% STRAT_classes defines the initial CryoGrid stratigraphy of stratigraphy 
% classes.
% Each stratigraphy class is identified by its name and a class_index 
% (integer number)  which makes it possible to define the same stratigraphy
% classs everal times with different parameters. For stratigraphy classes 
% interacting with a snow cover, a SNOW class must be % defined.
% Sleeping classes are defined within the stratigraphy are initialized, 
% but not included in the initial CryoGrid stratigraphy. Instead, they 
% are made accessible when requested by a trigger.
% An example is the formation of a lake on initially dry ground. A LAKE
% class is initialized as sleeping class, and added on top of the
% stratigraphy when sufficient water has accumulated.
% S. Westermann, T. Ingeman-Nielsen, J. Scheer, October 2020
%========================================================================

classdef STRAT_classes_variable_upper_layer < matlab.mixin.Copyable

    
    properties
		strat_classes_index
		PARA
        CONST
        depth
        class_name
        class_index
        snow_class
        sleeping_classes
    end
    
    methods
        
		
		function self = provide_PARA(self)
            self.PARA.thickness_upper_layer = [];	
			self.PARA.classes = [];
			self.PARA.snow_class_name = [];
			self.PARA.snow_class_index = [];
            self.PARA.sleeping_classes_name = [];
            self.PARA.sleeping_classes_index = [];
        end
		
        function self = provide_CONST(self)

        end
        
        function self = provide_STATVAR(self)

        end       
        
        function stratigraphy = finalize_init(stratigraphy, tile)
            stratigraphy.PARA.classes.depth(2,1) = stratigraphy.PARA.thickness_upper_layer;
            i=1;
            while i <= size(stratigraphy.PARA.classes.depth,1)-1
                if  stratigraphy.PARA.classes.depth(i+1,1) <=  stratigraphy.PARA.classes.depth(i,1)+0.2 %remove  class if layerThickness <0.2m
                    stratigraphy.PARA.classes.depth(i,:) = [];
                    if i==1
                        stratigraphy.PARA.classes.depth(i,1) = 0;
                    end
                    stratigraphy.PARA.classes.class_name(i,:) = [];
                    stratigraphy.PARA.classes.class_index(i,:) = [];
                end
                i=i+1;
            end
        end 
        
        
        %-------------param file generation-----
         function stratigraphy = param_file_info(stratigraphy)
             stratigraphy = provide_PARA(stratigraphy);
             %default
             stratigraphy.PARA.default_value = [];
             
             stratigraphy.PARA.STATVAR = [];
             stratigraphy.PARA.class_category = 'STRATIGRAPHY_CLASSES';

             stratigraphy.PARA.options.classes.name = 'STRAT_MATRIX';
             stratigraphy.PARA.options.classes.entries_x = {'class_name', 'class_index'};
             stratigraphy.PARA.options.classes.entries_y = {0; 20};
             stratigraphy.PARA.comment.classes = {'stratigraphy of subsurface classes'};
             
             stratigraphy.PARA.comment.snow_class_name = {'snow class to be used in the model run'};
             
             stratigraphy.PARA.options.sleeping_classes_name.name = 'H_LIST';
             stratigraphy.PARA.options.sleeping_classes_name.entries_x = {};
             stratigraphy.PARA.comment.sleeping_classes_name = {'classes that are added to the stratigraphy during the model run'};
             
             stratigraphy.PARA.options.sleeping_classes_index.name = 'H_LIST';
             stratigraphy.PARA.options.sleeping_classes_index.entries_x = {};
         end
        
        
        
        

    end
    
end

