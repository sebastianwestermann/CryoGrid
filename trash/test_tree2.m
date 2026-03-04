
CG_tree = append_next_level(provider, class(provider.RUN_INFO_CLASS), provider.RUN_INFO_CLASS.PARA.class_index, 1);
fig = visualizeTreeStruct(provider, CG_tree, 'Title', 'Simulation Parameters');


function res = append_next_level(provider, class_name, class_index, depth_index)

    res = [];
    fn2 = fieldnames(provider.CLASSES);
    index = find(strcmp(fn2, class_name));
    if isempty(index)
        disp(['class ' class_name  ' with index ' class_index ' is not defined'])
    else
        selected_class = provider.CLASSES.(fn2{index,1}){class_index,1};
        fn = fieldnames(selected_class.PARA);
        pos_of_classes=[];
        for i=1:size(fn,1)
            if size(fn{i,1},2)> 6 && strcmp(fn{i,1}(end-5:end), '_index') && ~strcmp(fn{i,1}, 'class_index')
                index_of_class = find(strcmp(fn, fn{i,1}(1:end-6)) | strcmp(fn, [fn{i,1}(1:end-6) '_name']));
                pos_of_classes = [pos_of_classes; [i index_of_class]];
                if ~isempty(selected_class.PARA.(fn{index_of_class,1}))
                    if ~iscell(selected_class.PARA.(fn{index_of_class,1}))
                        if ~isnan(selected_class.PARA.(fn{index_of_class,1}))
                            %recursive call to function
                            a = append_next_level(provider, selected_class.PARA.(fn{index_of_class,1}), selected_class.PARA.(fn{i,1}), depth_index+1);
                            res.(class(selected_class)).([fn{index_of_class,1} '__' selected_class.PARA.(fn{index_of_class,1})]) = a.(selected_class.PARA.(fn{index_of_class,1}));
                            % res.(class(selected_class)).(fn{index_of_class,1}) = ...
                            %     append_next_level(provider, selected_class.PARA.(fn{index_of_class,1}), selected_class.PARA.(fn{i,1}), depth_index+1);
                        else
                            res.(class(selected_class)).(fn{index_of_class,1}) = selected_class.PARA.(fn{index_of_class,1});
                        end
                    else %cell array
                        for j=1:size(selected_class.PARA.(fn{index_of_class,1}),1)
                            if ~isnan(selected_class.PARA.(fn{index_of_class,1}){j,1})
                                a = append_next_level(provider, selected_class.PARA.(fn{index_of_class,1}){j,1}, selected_class.PARA.(fn{i,1})(j,1), depth_index+1);
                                res.(class(selected_class)).([fn{index_of_class,1} '_' num2str(j) '__' selected_class.PARA.(fn{index_of_class,1}){j,1}]) = a.(selected_class.PARA.(fn{index_of_class,1}){j,1});
                            else
                                res.(class(selected_class)).([fn{index_of_class,1} '_' num2str(j) '__' selected_class.PARA.(fn{index_of_class,1}){j,1}]) = selected_class.PARA.(fn{index_of_class,1}){j,1};
                            end
                        end
                    end                    
                else
                    res.(class(selected_class)).(fn{index_of_class,1}) = [];
                end
               

            end
        end

        for i=1:size(fn,1)
            if ~any(i==pos_of_classes(:))
                if ~iscell(selected_class.PARA.(fn{i,1}))
                    if isstruct(selected_class.PARA.(fn{i,1})) %check whether there are classes in a struct
                        fn_struct = fieldnames(selected_class.PARA.(fn{i,1}));
                        index_of_class2 = [];
                        pos_of_classes2 = [];
                        for j=1:size(fn_struct,1)
                            if size(fn_struct{j,1},2)>6 && strcmp(fn_struct{j,1}(end-5:end), '_index')
                                index_of_class2 = find(strcmp(fn_struct, fn_struct{j,1}(1:end-6)) | strcmp(fn_struct, [fn_struct{j,1}(1:end-6) '_name']));
                                pos_of_classes2 = [pos_of_classes2; [j index_of_class2]];
                                if ~isempty(selected_class.PARA.(fn{i,1}).(fn_struct{index_of_class2,1}))
                                    if ~iscell(selected_class.PARA.(fn{i,1}).(fn_struct{index_of_class2,1}))
                                        if ~isnan(selected_class.PARA.(fn{i,1}).(fn_struct{index_of_class2,1}))
                                            %recursive call to function
                                            a = append_next_level(provider, selected_class.PARA.(fn{i,1}).(fn_struct{index_of_class2,1}), selected_class.PARA.(fn{i,1}).(fn_struct{j,1}), depth_index+1);
                                            res.(class(selected_class)).(fn{i,1}).(fn_struct{index_of_class2,1}) = a.(selected_class.PARA.(fn{i,1}).(fn_struct{index_of_class2,1}));
                                        else
                                            res.(class(selected_class)).(fn{i,1}).(fn_struct{index_of_class2,1}) = selected_class.PARA.(fn{i,1}).(fn_struct{index_of_class2,1});
                                        end
                                    else %cell array
                                        if size(selected_class.PARA.(fn{i,1}).(fn_struct{index_of_class2,1}),1)==0
                                            res.(class(selected_class)).(fn{i,1}).(fn_struct{index_of_class2,1}) = [];
                                        end
                                        for k=1:size(selected_class.PARA.(fn{i,1}).(fn_struct{index_of_class2,1}),1)
                                            if ~isnan(selected_class.PARA.(fn{i,1}).(fn_struct{index_of_class2,1}){k,1})
                                                a = append_next_level(provider, selected_class.PARA.(fn{i,1}).(fn_struct{index_of_class2,1}){k,1}, selected_class.PARA.(fn{i,1}).(fn_struct{j,1})(k,1), depth_index+1);
                                                res.(class(selected_class)).(fn{i,1}).([fn_struct{index_of_class2,1} '_' num2str(k) '_' selected_class.PARA.(fn{i,1}).(fn_struct{index_of_class2,1}){k,1}]) = a.(selected_class.PARA.(fn{i,1}).(fn_struct{index_of_class2,1}){k,1});
                                            else
                                                res.(class(selected_class)).(fn{i,1}).([fn_struct{index_of_class2,1} '_' num2str(k) '_' selected_class.PARA.(fn{i,1}).(fn_struct{index_of_class2,1}){k,1}]) = selected_class.PARA.(fn{i,1}).(fn_struct{index_of_class2,1}){k,1};
                                            end
                                        end
                                    end
                                else
                                    res.(class(selected_class)).(fn{i,1}).(fn_struct{index_of_class2,1}) = selected_class.PARA.(fn{i,1}).(fn_struct{index_of_class2,1});
                                end
                            end
                        end
                        for j=1:size(fn_struct,1)
                            if ~any(j==pos_of_classes2(:))
                                res.(class(selected_class)).(fn{i,1}).(fn_struct{j,1}) = selected_class.PARA.(fn{i,1}).(fn_struct{j,1});
                            end
                        end
                        %add go through the non-classes paras
                    else
                        res.(class(selected_class)).(fn{i,1}) = selected_class.PARA.(fn{i,1});
                    end
                else
                    for j=1:size(selected_class.PARA.(fn{i,1}),1)
                        res.(class(selected_class)).([fn{i,1} '_' num2str(j)]) = selected_class.PARA.(fn{i,1}){j,1};
                    end
                end
            end
        end
    end
end


%%
params = provider.CLASSES;
list_of_classes = [];
fn = fieldnames(params);
for i=1:size(fn,1)
    current_value = params.(fn{i,1});
    if size(current_value,1)==1
        fn2=fieldnames(current_value{1,1}.PARA);
        for k=1:size(fn2,1)
            if iscell(current_value{1,1}.PARA.(fn2{k,1}))
                for l=1:size(current_value{1,1}.PARA.(fn2{k,1}),1)
                    list_of_classes.(class(current_value{1,1})).(fn2{k,1}).(['entry_' num2str(l)]) = current_value{1,1}.PARA.(fn2{k,1}){l,1};
                end
            else
                list_of_classes.(class(current_value{1,1})).(fn2{k,1}) = current_value{1,1}.PARA.(fn2{k,1});
            end
        end
    else
        for j=1:size(current_value,1)
            %list_of_classes.(class(current_value{1,1})).([class(current_value{j,1}) '_' num2str(j)]) = current_value{j,1}.PARA;

            fn2=fieldnames(current_value{j,1}.PARA);
            for k=1:size(fn2,1)
                if iscell(current_value{j,1}.PARA.(fn2{k,1}))
                    for l=1:size(current_value{j,1}.PARA.(fn2{k,1}),1)
                        list_of_classes.(class(current_value{1,1})).([class(current_value{j,1}) '_' num2str(j)]).(fn2{k,1}).(['entry_' num2str(l)]) = ...
                            current_value{j,1}.PARA.(fn2{k,1}){l,1};
                    end
                else
                    list_of_classes.(class(current_value{1,1})).([class(current_value{j,1}) '_' num2str(j)]).(fn2{k,1}) = current_value{j,1}.PARA.(fn2{k,1});
                end
            end

        end
    end
end


%%

% Launch visualizer
fig = visualizeTreeStruct(provider, list_of_classes, 'Title', 'Simulation Parameters');

% Later, retrieve modified parameters:
% modifiedParams = fig.UserData.data;


function fig = visualizeTreeStruct(provider, treeStruct, varargin)
%VISUALIZETREESTRUCT Interactive visualization and editor for nested structs

p = inputParser;
addParameter(p, 'ReadOnly', false, @islogical);
addParameter(p, 'Title', 'Tree Structure Editor', @ischar);
parse(p, varargin{:});

opts = p.Results;

% Create main figure
fig = uifigure('Name', opts.Title, 'Position', [100 100 900 600]);
fig.UserData = struct('data', treeStruct, 'modified', false);

% Main layout: 2 rows (content + buttons), 1 column
mainLayout = uigridlayout(fig, [2 1]);
mainLayout.RowHeight = {'1x', 50};
mainLayout.Padding = [5 5 5 5];

% Content grid: Tree (left 1/3) | Properties (right 2/3)
contentGrid = uigridlayout(mainLayout, [1 2]);
contentGrid.Layout.Row = 1;
contentGrid.Layout.Column = 1;
contentGrid.ColumnWidth = {'1x', '2x'};
contentGrid.Padding = [0 0 0 0];

%% Left Panel: Tree View
treePanel = uipanel(contentGrid, 'Title', 'Hierarchy');
treePanel.Layout.Row = 1;
treePanel.Layout.Column = 1;

treeGrid = uigridlayout(treePanel, [1 1]);
treeGrid.Padding = [0 0 0 0];

tree = uitree(treeGrid, 'SelectionChangedFcn', @(src, event) nodeSelected(src, event));
tree.Layout.Row = 1;
tree.Layout.Column = 1;

% Build tree recursively
rootName = inputname(1);
if isempty(rootName), rootName = 'Root'; end
root = uitreenode(tree, 'Text', rootName, ...
    'NodeData', struct('path', {{}}, 'isStruct', true));

buildTreeNodes(root, treeStruct, {});
expand(root);

%% Right Panel: Property Editor - TWO COLUMN LAYOUT
propPanel = uipanel(contentGrid, 'Title', 'Parameters');
propPanel.Layout.Row = 1;
propPanel.Layout.Column = 2;

propGrid = uigridlayout(propPanel);
propGrid.Padding = [10 10 10 10];

%% Bottom button bar
btnGrid = uigridlayout(mainLayout, [1 3]);
btnGrid.Layout.Row = 2;
btnGrid.Layout.Column = 1;
btnGrid.RowHeight = {'1x'};
btnGrid.ColumnWidth = {'1x', 100, 100};
btnGrid.Padding = [5 5 5 5];

statusLabel = uilabel(btnGrid, 'Text', 'Select a node to view/edit parameters');
statusLabel.Layout.Row = 1;
statusLabel.Layout.Column = 1;

exportBtn = uibutton(btnGrid, 'push', 'Text', 'Export', ...
    'ButtonPushedFcn', @(~,~) exportData(fig));
exportBtn.Layout.Row = 1;
exportBtn.Layout.Column = 2;

resetBtn = uibutton(btnGrid, 'push', 'Text', 'Reset', ...
    'ButtonPushedFcn', @(~,~) resetData(fig, treeStruct));
resetBtn.Layout.Row = 1;
resetBtn.Layout.Column = 3;

if opts.ReadOnly
    exportBtn.Enable = 'off';
    resetBtn.Enable = 'off';
end

    %% Nested Functions
    
    function buildTreeNodes(parent, s, path)
        fields = fieldnames(s);
        for i = 1:numel(fields)
            fname = fields{i};
            val = s.(fname);
            currentPath = [path {fname}];
            
            if isstruct(val)
                txt = sprintf('%s \x25B6', fname); 
                node = uitreenode(parent, 'Text', txt, ...
                    'NodeData', struct('path', {currentPath}, 'isStruct', true));
                buildTreeNodes(node, val, currentPath);
            else
                uitreenode(parent, 'Text', fname, ...
                    'NodeData', struct('path', {currentPath}, 'isStruct', false));
            end
        end
    end
    
    function val = getValueFromPath(path)
        val = fig.UserData.data;
        for i = 1:numel(path)
            val = val.(path{i});
        end
    end
    
    function displayNode(data)
        delete(propGrid.Children);
        
        if data.isStruct
            createStructEditor(data);
        else
            createLeafEditor(data);
        end
    end
    
    function nodeSelected(~, event)
        if isempty(event.SelectedNodes)
            return;
        end
        
        nodeInfo = event.SelectedNodes(1).NodeData;
        currentValue = getValueFromPath(nodeInfo.path);
        
        displayData = struct();
        displayData.path = nodeInfo.path;
        displayData.value = currentValue;
        displayData.isStruct = nodeInfo.isStruct;
        
        if ~nodeInfo.isStruct
            displayData.field = nodeInfo.path{end};
        end
        
        displayNode(displayData);
    end
    
    function createStructEditor(data)
        fields = fieldnames(data.value);
        n = numel(fields);
        
        if n == 0
            uilabel(propGrid, 'Text', '(Empty struct)', ...
                'HorizontalAlignment', 'center');
            propGrid.RowHeight = {'1x'};
            propGrid.ColumnWidth = {'1x'};
            return;
        end
        
        % Calculate rows needed (2 parameters per row)
        numRows = ceil(n / 2);
        
        % Grid layout: Label | Value | Type | Label | Value | Type
        % Values now 100px (double width), labels stay flexible ('1x'), types 40px
        propGrid.RowHeight = repmat({35}, 1, numRows);
        propGrid.ColumnWidth = {'1x', 100, 40, '1x', 100, 40};
        
        for i = 1:n
            % Calculate position in 2-column layout
            row = ceil(i / 2);
            isLeft = (mod(i, 2) == 1);
            colOffset = 0;
            if ~isLeft
                colOffset = 3;
            end
            
            fname = fields{i};
            val = data.value.(fname);
            
            % Label - flexible width
            lbl = uilabel(propGrid, 'Text', fname, ...
                'HorizontalAlignment', 'right', ...
                'FontWeight', 'bold');
            lbl.Layout.Row = row;
            lbl.Layout.Column = 1 + colOffset;
            
            % Value/Editor
            if isstruct(val)
                btnPath = [data.path {fname}];
                btn = uibutton(propGrid, 'push', 'Text', 'Struct...', ...
                    'ButtonPushedFcn', @(~,~) selectNodeByPath(btnPath));
                btn.Layout.Row = row;
                btn.Layout.Column = 2 + colOffset;
                
                typeLbl = uilabel(propGrid, 'Text', '(struct)');
                typeLbl.Layout.Row = row;
                typeLbl.Layout.Column = 3 + colOffset;
            else
                createEditorControl(row, 2 + colOffset, data.path, fname, val);
                
                typeLbl = uilabel(propGrid, 'Text', shortType(val));
                typeLbl.Layout.Row = row;
                typeLbl.Layout.Column = 3 + colOffset;
            end
        end
    end
    
    function createLeafEditor(data)
        % Use same layout with wider values
        propGrid.RowHeight = {35};
        propGrid.ColumnWidth = {'1x', 100, 40, '1x', 100, 40};
        
        lbl = uilabel(propGrid, 'Text', data.field, ...
            'HorizontalAlignment', 'right', ...
            'FontWeight', 'bold');
        lbl.Layout.Row = 1;
        lbl.Layout.Column = 1;
        
        if numel(data.path) > 1
            pathToParent = data.path(1:end-1);
        else
            pathToParent = {};
        end
        fieldName = data.path{end};
        
        createEditorControl(1, 2, pathToParent, fieldName, data.value);
        
        typeLbl = uilabel(propGrid, 'Text', shortType(data.value));
        typeLbl.Layout.Row = 1;
        typeLbl.Layout.Column = 3;
    end
    
    function createEditorControl(row, col, path, field, val)
        created = false;
        
        % Check if value is editable numeric (scalar, finite, not NaN)
        isEditableNumeric = isnumeric(val) && isscalar(val) && ~isnan(val) && ~isinf(val);
        
        if ~opts.ReadOnly && isEditableNumeric
            created = true;
            if isinteger(val)
                edit = uispinner(propGrid, 'Value', double(val), ...
                    'Limits', [-inf inf], ...
                    'ValueChangedFcn', @(src,~) updateStruct(path, field, cast(src.Value, class(val))));
            else
                edit = uieditfield(propGrid, 'numeric', 'Value', val, ...
                    'ValueChangedFcn', @(src,~) updateStruct(path, field, src.Value));
            end
            edit.Layout.Row = row;
            edit.Layout.Column = col;
            
        elseif ~opts.ReadOnly && (ischar(val) || isstring(val))
            created = true;
            strVal = char(val);
            edit = uieditfield(propGrid, 'text', 'Value', strVal, ...
                'ValueChangedFcn', @(src,~) updateStruct(path, field, char(src.Value)));
            edit.Layout.Row = row;
            edit.Layout.Column = col;
            
        elseif ~opts.ReadOnly && islogical(val) && isscalar(val)
            created = true;
            edit = uicheckbox(propGrid, 'Text', '', 'Value', val, ...
                'ValueChangedFcn', @(src,~) updateStruct(path, field, src.Value));
            edit.Layout.Row = row;
            edit.Layout.Column = col;
        end
        
        if created
            if opts.ReadOnly, edit.Editable = false; end
        else
            % Display as read-only text for non-editable values
            txt = uilabel(propGrid, 'Text', valueSummary(val));
            txt.HorizontalAlignment = 'left';
            txt.Layout.Row = row;
            txt.Layout.Column = col;
        end
    end
    
    function str = shortType(v)
        sz = size(v);
        if isscalar(v)
            str = class(v);
        else
            str = sprintf('[%s]', strjoin(string(sz), '×'));
        end
    end
    
    function str = valueSummary(v)
        if isnumeric(v)
            if isempty(v)
                str = '[]';
            elseif isscalar(v)
                if isnan(v)
                    str = 'NaN';
                elseif isinf(v)
                    str = 'Inf';
                else
                    str = num2str(v, 4);
                end
            else
                str = sprintf('[%s]', mat2str(size(v)));
            end
        elseif ischar(v) || isstring(v)
            str = sprintf('"%s"', char(v));
        elseif islogical(v)
            str = mat2str(v);
        else
            str = sprintf('<%s>', class(v));
        end
        if length(str) > 20
            str = [str(1:17) '...'];
        end
    end
    
    function selectNodeByPath(targetPath)
        node = findNodeInTree(root, targetPath);
        if ~isempty(node)
            tree.SelectedNodes = node;
            expand(node);
            nodeInfo = node.NodeData;
            currentValue = getValueFromPath(nodeInfo.path);
            displayData = struct();
            displayData.path = nodeInfo.path;
            displayData.value = currentValue;
            displayData.isStruct = nodeInfo.isStruct;
            if ~nodeInfo.isStruct
                displayData.field = nodeInfo.path{end};
            end
            displayNode(displayData);
        end
    end
    
    function node = findNodeInTree(parent, targetPath)
        node = [];
        if isempty(parent.Children)
            return;
        end
        
        children = parent.Children;
        for i = 1:numel(children)
            childPath = children(i).NodeData.path;
            if isequal(childPath, targetPath)
                node = children(i);
                return;
            elseif numel(childPath) < numel(targetPath) && ...
                   isequal(childPath, targetPath(1:numel(childPath)))
                node = findNodeInTree(children(i), targetPath);
                if ~isempty(node), return; end
            end
        end
    end
    
    function updateStruct(path, field, newVal)
        data = fig.UserData.data;
        data = setfield(data, path{:}, field, newVal);
        fig.UserData.data = data;
        fig.UserData.modified = true;
    end
    
    function updateArray(path, field, strVal)
        try
            newVal = eval(strVal);
            updateStruct(path, field, newVal);
            uialert(fig, 'Array updated successfully', 'Success', 'Icon', 'success', 'Modal', false);
        catch ME
            uialert(fig, sprintf('Invalid array syntax: %s', ME.message), 'Error', 'Icon', 'error', 'Modal', false);
        end
    end
    
    function exportData(~)
        assignin('base', 'treeParams', fig.UserData.data);
        uialert(fig, 'Exported to base workspace as ''treeParams''', 'Success', 'Icon', 'success', 'Modal', false);
    end

    function resetData(fig, original)
        fig.UserData.data = original;
        fig.UserData.modified = false;
        delete(root.Children);
        buildTreeNodes(root, fig.UserData.data, {});
        expand(root);
        tree.SelectedNodes = root;
        displayNode(struct('path', {{}}, 'value', fig.UserData.data, 'isStruct', true));
    end

%% Initialize: Select root and display its contents
tree.SelectedNodes = root;
displayNode(struct('path', {{}}, 'value', fig.UserData.data, 'isStruct', true));
end





% function fig = visualizeTreeStruct(provider, treeStruct, varargin)
% %VISUALIZETREESTRUCT Interactive visualization and editor for nested structs
% 
% p = inputParser;
% addParameter(p, 'ReadOnly', false, @islogical);
% addParameter(p, 'Title', 'Tree Structure Editor', @ischar);
% parse(p, varargin{:});
% 
% opts = p.Results;
% 
% % Create main figure
% fig = uifigure('Name', opts.Title, 'Position', [100 100 900 600]);
% fig.UserData = struct('data', treeStruct, 'modified', false);
% 
% % Main layout: 2 rows (content + buttons), 1 column
% mainLayout = uigridlayout(fig, [2 1]);
% mainLayout.RowHeight = {'1x', 50};
% mainLayout.Padding = [5 5 5 5];
% 
% % Content grid: Tree (left 1/3) | Properties (right 2/3)
% contentGrid = uigridlayout(mainLayout, [1 2]);
% contentGrid.Layout.Row = 1;
% contentGrid.Layout.Column = 1;
% contentGrid.ColumnWidth = {'1x', '2x'};
% contentGrid.Padding = [0 0 0 0];
% 
% %% Left Panel: Tree View
% treePanel = uipanel(contentGrid, 'Title', 'Hierarchy');
% treePanel.Layout.Row = 1;
% treePanel.Layout.Column = 1;
% 
% treeGrid = uigridlayout(treePanel, [1 1]);
% treeGrid.Padding = [0 0 0 0];
% 
% tree = uitree(treeGrid, 'SelectionChangedFcn', @(src, event) nodeSelected(src, event));
% tree.Layout.Row = 1;
% tree.Layout.Column = 1;
% 
% % Build tree recursively
% rootName = inputname(1);
% if isempty(rootName), rootName = 'Root'; end
% root = uitreenode(tree, 'Text', rootName, ...
%     'NodeData', struct('path', {{}}, 'isStruct', true));
% 
% buildTreeNodes(root, treeStruct, {});
% expand(root);
% 
% %% Right Panel: Property Editor - TWO COLUMN LAYOUT
% propPanel = uipanel(contentGrid, 'Title', 'Parameters');
% propPanel.Layout.Row = 1;
% propPanel.Layout.Column = 2;
% 
% propGrid = uigridlayout(propPanel);
% propGrid.Padding = [10 10 10 10];
% 
% %% Bottom button bar
% btnGrid = uigridlayout(mainLayout, [1 3]);
% btnGrid.Layout.Row = 2;
% btnGrid.Layout.Column = 1;
% btnGrid.RowHeight = {'1x'};
% btnGrid.ColumnWidth = {'1x', 100, 100};
% btnGrid.Padding = [5 5 5 5];
% 
% statusLabel = uilabel(btnGrid, 'Text', 'Select a node to view/edit parameters');
% statusLabel.Layout.Row = 1;
% statusLabel.Layout.Column = 1;
% 
% exportBtn = uibutton(btnGrid, 'push', 'Text', 'Export', ...
%     'ButtonPushedFcn', @(~,~) exportData(fig));
% exportBtn.Layout.Row = 1;
% exportBtn.Layout.Column = 2;
% 
% resetBtn = uibutton(btnGrid, 'push', 'Text', 'Reset', ...
%     'ButtonPushedFcn', @(~,~) resetData(fig, treeStruct));
% resetBtn.Layout.Row = 1;
% resetBtn.Layout.Column = 3;
% 
% if opts.ReadOnly
%     exportBtn.Enable = 'off';
%     resetBtn.Enable = 'off';
% end
% 
%     %% Nested Functions
% 
%     function buildTreeNodes(parent, s, path)
%         fields = fieldnames(s);
%         for i = 1:numel(fields)
%             fname = fields{i};
%             val = s.(fname);
%             currentPath = [path {fname}];
% 
%             if isstruct(val)
%                 txt = sprintf('%s \x25B6', fname); 
%                 node = uitreenode(parent, 'Text', txt, ...
%                     'NodeData', struct('path', {currentPath}, 'isStruct', true));
%                 buildTreeNodes(node, val, currentPath);
%             else
%                 uitreenode(parent, 'Text', fname, ...
%                     'NodeData', struct('path', {currentPath}, 'isStruct', false));
%             end
%         end
%     end
% 
%     function val = getValueFromPath(path)
%         val = fig.UserData.data;
%         for i = 1:numel(path)
%             val = val.(path{i});
%         end
%     end
% 
%     function displayNode(data)
%         delete(propGrid.Children);
% 
%         if data.isStruct
%             createStructEditor(data);
%         else
%             createLeafEditor(data);
%         end
%     end
% 
%     function nodeSelected(~, event)
%         if isempty(event.SelectedNodes)
%             return;
%         end
% 
%         nodeInfo = event.SelectedNodes(1).NodeData;
%         currentValue = getValueFromPath(nodeInfo.path);
% 
%         displayData = struct();
%         displayData.path = nodeInfo.path;
%         displayData.value = currentValue;
%         displayData.isStruct = nodeInfo.isStruct;
% 
%         if ~nodeInfo.isStruct
%             displayData.field = nodeInfo.path{end};
%         end
% 
%         displayNode(displayData);
%     end
% 
%     function createStructEditor(data)
%         fields = fieldnames(data.value);
%         n = numel(fields);
% 
%         if n == 0
%             uilabel(propGrid, 'Text', '(Empty struct)', ...
%                 'HorizontalAlignment', 'center');
%             propGrid.RowHeight = {'1x'};
%             propGrid.ColumnWidth = {'1x'};
%             return;
%         end
% 
%         % Calculate rows needed (2 parameters per row)
%         numRows = ceil(n / 2);
% 
%         % Grid layout: Label | Value | Type | Label | Value | Type
%         propGrid.RowHeight = repmat({35}, 1, numRows);
%         propGrid.ColumnWidth = {100, '1x', 70, 100, '1x', 70};
% 
%         for i = 1:n
%             % Calculate position in 2-column layout
%             row = ceil(i / 2);
%             isLeft = (mod(i, 2) == 1);
%             colOffset = 0;
%             if ~isLeft
%                 colOffset = 3;
%             end
% 
%             fname = fields{i};
%             val = data.value.(fname);
% 
%             % Label
%             lbl = uilabel(propGrid, 'Text', fname, ...
%                 'HorizontalAlignment', 'right', ...
%                 'FontWeight', 'bold');
%             lbl.Layout.Row = row;
%             lbl.Layout.Column = 1 + colOffset;
% 
%             % Value/Editor
%             if isstruct(val)
%                 btnPath = [data.path {fname}];
%                 btn = uibutton(propGrid, 'push', 'Text', 'Struct...', ...
%                     'ButtonPushedFcn', @(~,~) selectNodeByPath(btnPath));
%                 btn.Layout.Row = row;
%                 btn.Layout.Column = 2 + colOffset;
% 
%                 typeLbl = uilabel(propGrid, 'Text', '(struct)');
%                 typeLbl.Layout.Row = row;
%                 typeLbl.Layout.Column = 3 + colOffset;
%             else
%                 createEditorControl(row, 2 + colOffset, data.path, fname, val);
% 
%                 typeLbl = uilabel(propGrid, 'Text', shortType(val));
%                 typeLbl.Layout.Row = row;
%                 typeLbl.Layout.Column = 3 + colOffset;
%             end
%         end
%     end
% 
%     function createLeafEditor(data)
%         propGrid.RowHeight = {35};
%         propGrid.ColumnWidth = {100, '1x', 70, 100, '1x', 70};
% 
%         lbl = uilabel(propGrid, 'Text', data.field, ...
%             'HorizontalAlignment', 'right', ...
%             'FontWeight', 'bold');
%         lbl.Layout.Row = 1;
%         lbl.Layout.Column = 1;
% 
%         if numel(data.path) > 1
%             pathToParent = data.path(1:end-1);
%         else
%             pathToParent = {};
%         end
%         fieldName = data.path{end};
% 
%         createEditorControl(1, 2, pathToParent, fieldName, data.value);
% 
%         typeLbl = uilabel(propGrid, 'Text', shortType(data.value));
%         typeLbl.Layout.Row = 1;
%         %typeLbl.Layout = 3;
%         typeLbl.Layout.Column = 3;
%     end
% 
%     function createEditorControl(row, col, path, field, val)
%         created = false;
% 
%         % Check if value is editable numeric (scalar, finite, not NaN)
%         isEditableNumeric = isnumeric(val) && isscalar(val) && ~isnan(val) && ~isinf(val);
% 
%         if ~opts.ReadOnly && isEditableNumeric
%             created = true;
%             if isinteger(val)
%                 edit = uispinner(propGrid, 'Value', double(val), ...
%                     'Limits', [-inf inf], ...
%                     'ValueChangedFcn', @(src,~) updateStruct(path, field, cast(src.Value, class(val))));
%             else
%                 edit = uieditfield(propGrid, 'numeric', 'Value', val, ...
%                     'ValueChangedFcn', @(src,~) updateStruct(path, field, src.Value));
%             end
%             edit.Layout.Row = row;
%             edit.Layout.Column = col;
% 
%         elseif ~opts.ReadOnly && (ischar(val) || isstring(val))
%             created = true;
%             strVal = char(val);
%             edit = uieditfield(propGrid, 'text', 'Value', strVal, ...
%                 'ValueChangedFcn', @(src,~) updateStruct(path, field, char(src.Value)));
%             edit.Layout.Row = row;
%             edit.Layout.Column = col;
% 
%         elseif ~opts.ReadOnly && islogical(val) && isscalar(val)
%             created = true;
%             edit = uicheckbox(propGrid, 'Text', '', 'Value', val, ...
%                 'ValueChangedFcn', @(src,~) updateStruct(path, field, src.Value));
%             edit.Layout.Row = row;
%             edit.Layout.Column = col;
%         end
% 
%         if created
%             if opts.ReadOnly, edit.Editable = false; end
%         else
%             % Display as read-only text for non-editable values
%             txt = uilabel(propGrid, 'Text', valueSummary(val));
%             txt.HorizontalAlignment = 'left';
%             txt.Layout.Row = row;
%             txt.Layout.Column = col;
%         end
%     end
% 
%     function str = shortType(v)
%         sz = size(v);
%         if isscalar(v)
%             str = class(v);
%         else
%             str = sprintf('[%s]', strjoin(string(sz), '×'));
%         end
%     end
% 
%     function str = valueSummary(v)
%         if isnumeric(v)
%             if isempty(v)
%                 str = '[]';
%             elseif isscalar(v)
%                 if isnan(v)
%                     str = 'NaN';
%                 elseif isinf(v)
%                     str = 'Inf';
%                 else
%                     str = num2str(v, 4);
%                 end
%             else
%                 str = sprintf('[%s]', mat2str(size(v)));
%             end
%         elseif ischar(v) || isstring(v)
%             str = sprintf('"%s"', char(v));
%         elseif islogical(v)
%             str = mat2str(v);
%         else
%             str = sprintf('<%s>', class(v));
%         end
%         if length(str) > 20
%             str = [str(1:17) '...'];
%         end
%     end
% 
%     function selectNodeByPath(targetPath)
%         node = findNodeInTree(root, targetPath);
%         if ~isempty(node)
%             tree.SelectedNodes = node;
%             expand(node);
%             nodeInfo = node.NodeData;
%             currentValue = getValueFromPath(nodeInfo.path);
%             displayData = struct();
%             displayData.path = nodeInfo.path;
%             displayData.value = currentValue;
%             displayData.isStruct = nodeInfo.isStruct;
%             if ~nodeInfo.isStruct
%                 displayData.field = nodeInfo.path{end};
%             end
%             displayNode(displayData);
%         end
%     end
% 
%     function node = findNodeInTree(parent, targetPath)
%         node = [];
%         if isempty(parent.Children)
%             return;
%         end
% 
%         children = parent.Children;
%         for i = 1:numel(children)
%             childPath = children(i).NodeData.path;
%             if isequal(childPath, targetPath)
%                 node = children(i);
%                 return;
%             elseif numel(childPath) < numel(targetPath) && ...
%                    isequal(childPath, targetPath(1:numel(childPath)))
%                 node = findNodeInTree(children(i), targetPath);
%                 if ~isempty(node), return; end
%             end
%         end
%     end
% 
%     function updateStruct(path, field, newVal)
%         data = fig.UserData.data;
%         data = setfield(data, path{:}, field, newVal);
%         fig.UserData.data = data;
%         fig.UserData.modified = true;
%     end
% 
%     function updateArray(path, field, strVal)
%         try
%             newVal = eval(strVal);
%             updateStruct(path, field, newVal);
%             uialert(fig, 'Array updated successfully', 'Success', 'Icon', 'success', 'Modal', false);
%         catch ME
%             uialert(fig, sprintf('Invalid array syntax: %s', ME.message), 'Error', 'Icon', 'error', 'Modal', false);
%         end
%     end
% 
%     function exportData(~)
%         assignin('base', 'treeParams', fig.UserData.data);
%         uialert(fig, 'Exported to base workspace as ''treeParams''', 'Success', 'Icon', 'success', 'Modal', false);
%     end
% 
%     function resetData(fig, original)
%         fig.UserData.data = original;
%         fig.UserData.modified = false;
%         delete(root.Children);
%         buildTreeNodes(root, fig.UserData.data, {});
%         expand(root);
%         tree.SelectedNodes = root;
%         displayNode(struct('path', {{}}, 'value', fig.UserData.data, 'isStruct', true));
%     end
% 
% %% Initialize: Select root and display its contents
% tree.SelectedNodes = root;
% displayNode(struct('path', {{}}, 'value', fig.UserData.data, 'isStruct', true));
% end

