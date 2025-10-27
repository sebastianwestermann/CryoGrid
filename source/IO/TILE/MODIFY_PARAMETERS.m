classdef MODIFY_PARAMETERS < matlab.mixin.Copyable


    properties
        PARA
        CONST
        STATVAR
        TEMP
    end

    methods

        function mod = provide_PARA(mod)

            mod.PARA.modify_class_name = [];
            mod.PARA.modify_class_index = [];
            mod.PARA.variable_in_class = [];
            mod.PARA.index_V = [];
            mod.PARA.index_H = [];
            mod.PARA.updated_value = [];

        end

        function mod = provide_CONST(mod)

        end

        function mod = provide_STATVAR(mod)

        end


        function mod = finalize_init(mod, tile)
            if isempty(mod.PARA.index_V) || (size(mod.PARA.index_V,1)==1 && isnan(mod.PARA.index_V))
                mod.PARA.index_V = ones(size(mod.PARA.modify_class_index,1));
            end
            mod.PARA.index_V(isnan(mod.PARA.index_V)) = 1;
            if isempty(mod.PARA.index_H) || (size(mod.PARA.index_H,1)==1 && isnan(mod.PARA.index_H))
                mod.PARA.index_H = ones(size(mod.PARA.modify_class_index,1));
            end
            mod.PARA.index_H(isnan(mod.PARA.index_H)) = 1;

        end



        function tile = modify(mod, tile)

            %establish pointers to all classes in the stratigraphy
            mod = set_pointers2classes(mod, tile);

            %rewrite PARA in PROVIDER and all classes in the
            %stratigraphy
            for i=1:size(mod.PARA.modify_class_name,1)
                tile.RUN_INFO.PPROVIDER.CLASSES.(mod.PARA.modify_class_name{i,1}){mod.PARA.modify_class_index(i,1),1}.PARA.(mod.PARA.variable_in_class{i,1})(mod.PARA.index_V(i,1), mod.PARA.index_H(i,1)) = ...
                    mod.PARA.updated_value(i,1);
                for j=1:size(mod.PARA.modify_class_pointer{i,1},1)
                    mod.PARA.modify_class_pointer{i,1}(j,1).PARA.(mod.PARA.variable_in_class{i,1})(mod.PARA.index_V(i,1), mod.PARA.index_H(i,1)) = mod.PARA.updated_value(i,1);
                end
            end
            mod.PARA.modify_class_pointer = [];

        end

        function mod = set_pointers2classes(mod, tile)

            for i=1:size(mod.PARA.modify_class_name,1)
                mod.PARA.modify_class_pointer{i,1} = {};
            end
            mod = find_classes_in_variable(mod, tile, 1);

        end

        function mod = find_classes_in_variable(mod, current_class, level)

            class_name = class(current_class);
            if isobject(current_class)
                for i=1:size(mod.PARA.modify_class_name, 1)
                    if strcmp(mod.PARA.modify_class_name{i,1}, class_name)
                        if  mod.PARA.modify_class_index(i,1) == current_class.PARA.class_index
                            dec = 1;
                            for j=1:size(mod.PARA.modify_class_pointer{i,1}, 1)
                                if isequal(mod.PARA.modify_class_pointer{i,1}(j,1), current_class)
                                    dec = 0;
                                end
                            end
                            if dec
                                mod.PARA.modify_class_pointer{i,1} = [mod.PARA.modify_class_pointer{i,1}; current_class];
                            end
                        end
                    end
                end
            end
            if isstruct(current_class) || isobject(current_class)
                variables = fieldnames(current_class);
                for i = 1:size(variables,1)
                    if iscell(current_class.(variables{i,1}))
                        if level <= 5
                            for j=1:size(current_class.(variables{i,1}), 1)
                                for k=1:size(current_class.(variables{i,1}), 2)
                                    mod = find_classes_in_variable(mod, current_class.(variables{i,1}){j,k}, level+1);
                                end
                            end
                        end
                    elseif isstruct(current_class.(variables{i,1}))
                        if level <= 5
                            mod = find_classes_in_variable(mod, current_class.(variables{i,1}), level+1);
                        end
                    elseif isobject(current_class.(variables{i,1})) && ~strcmp(variables{i,1}, 'RUN_INFO')
                        if level <= 5
                            mod = find_classes_in_variable(mod, current_class.(variables{i,1}), level+1);
                        end
                    end
                end
            end
        end

    end
end

