classdef MODIFY_CLASSES < matlab.mixin.Copyable


    properties
        PARA
        CONST
        STATVAR
        TEMP
    end

    methods

        function mod = provide_PARA(mod)

            mod.PARA.existing_class = [];
            mod.PARA.existing_class_index = [];
            mod.PARA.new_class = [];
            mod.PARA.new_class_index = [];

        end

        function mod = provide_CONST(mod)

        end

        function mod = provide_STATVAR(mod)

        end


        function mod = finalize_init(mod, tile)

        end



        function tile = modify(mod, tile)

            for i=1:size(mod.PARA.existing_class,1)
                new_class = copy(tile.RUN_INFO.PPROVIDER.CLASSES.(mod.PARA.new_class{i,1}){mod.PARA.new_class_index(i,1),1});
                ia_class = get_IA_class(mod.PARA.existing_class(i,1),  mod.PARA.new_class(i,1));

                CURRENT = tile.TOP.NEXT;
                while ~strcmp(class(CURRENT), mod.PARA.existing_class(i,1)) &&  ~strcmp(class(CURRENT), 'Bottom') %check also class number
                    CURRENT = CURRENT.NEXT;
                end
                if ~strcmp(class(CURRENT), 'Bottom')
                    CURRENT = change_class(ia_class, CURRENT, new_class, tile);
                end
            end

        end

    end
end

