%========================================================================
% CryoGrid INTERACTION (IA) for BGC classes to remove the FLIPFLOP
%========================================================================

classdef IA_FLIPFLOP2NORMAL_BGC <  IA_BASE
    
    properties

    end
    
    methods
        
        function new_class = change_class(ia, existing_class, new_class, tile)

            if strcmp(class(existing_class), 'GROUND_store_flip_flop_singleClass_BGC')
                existing_class.TEMP.next_flip_time = tile.t;
                existing_class = switch2model_flip_flop_BGC(existing_class, tile);
            end
            variables = fieldnames(existing_class);
            variables_new = fieldnames(new_class);
            for i=1:size(variables,1)
                if any(strcmp(variables_new, variables{i,1}))
                    new_class.(variables{i,1}) = existing_class.(variables{i,1});
                end
            end

            new_class.BGC = existing_class.BGC;
            new_class.BGC.PARENT = new_class;
            new_class.IA_BGC =  existing_class.IA_BGC;
            new_class.IA_BGC.BGC = new_class.BGC;
            new_class.IA_BGC.GROUND = new_class;
            new_class.BGC.IA_BGC = new_class.IA_BGC; 
                        
            new_class.PREVIOUS.NEXT = new_class;
            new_class.NEXT.PREVIOUS = new_class;

            ia_class = get_IA_class( class(new_class.PREVIOUS), class(new_class));
            new_class.IA_PREVIOUS = ia_class;
            new_class.PREVIOUS.IA_NEXT = ia_class;
            ia_class.NEXT = new_class;
            ia_class.PREVIOUS = new_class.PREVIOUS;
            finalize_init(new_class.IA_PREVIOUS, tile);

            ia_class = get_IA_class(class(new_class), class(new_class.NEXT));
            new_class.IA_NEXT = ia_class;
            new_class.NEXT.IA_PREVIOUS = ia_class;
            ia_class.PREVIOUS = new_class;
            ia_class.NEXT = new_class.NEXT;
            finalize_init(new_class.IA_NEXT, tile);

            existing_class.NEXT=[];
            existing_class.PREVIOUS=[];
            existing_class.IA_PREVIOUS=[];
            existing_class.IA_NEXT=[];

        end

    end
end