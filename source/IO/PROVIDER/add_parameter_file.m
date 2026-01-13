classdef add_parameter_file 
    
    properties
        PARA
        CONST
        STATVAR
    end

    methods

        function add_file = provide_PARA(add_file)
            add_file.PARA.parameter_file_folder = [];
            add_file.PARA.parameter_file_name = [];
        end

        function add_file = provide_CONST(add_file)

        end

        function add_file = provide_STATVAR(add_file)

        end

    end
end

