% 
% Below is a self-contained MATLAB function that visualizes a nested-struct (tree-like) and lets the user view and edit parameters at any level. It uses MATLAB's uitree/uitreenode for the tree view and a dynamic property panel on the right that shows editable controls for the selected node's non-struct fields. Changes are written back into the struct. The function returns the updated struct when the GUI is closed.
% 
% Usage:
% 
% Call updated = structEditorGui(origStruct). A figure opens. Click nodes to inspect/edit parameters. Close the figure to return the updated struct.
% Paste the function into a file named structEditorGui.m.

function updated = structEditorGui(s)
% structEditorGui  Visualize and edit a nested struct as a tree.
%
% updated = structEditorGui(s)
%
% Input:
%   s - nested struct (fields may be scalars, vectors, strings, numeric, logical,
%       cellstr etc. or another struct to represent nested levels)
%
% Output:
%   updated - modified struct returned after figure is closed.
%
% Features:
% - Left: tree view of the nested struct
% - Right: when selecting a node, editable controls are generated for all
%   non-struct fields at that node
% - Edits are propagated to the underlying struct
% - Use "Refresh tree" to rebuild tree if you changed structure programmatically
%
% Note: Requires MATLAB support for uitree/uitreenode (R2016b+ recommended).

if nargin==0
    s = struct('A', struct('a1', 1, 'a2', 'hello'), ...
               'B', struct('b1', true, 'b2', [1 2 3]), ...
               'C', 42);
end

% Use a handle to store current struct so nested functions can modify it
data.struct = s;

% Create figure
fig = uifigure('Name','Struct Editor','NumberTitle','off', ...
    'MenuBar','none','ToolBar','none', 'Position',[200 200 800 500], ...
    'CloseRequestFcn',@onClose);

% Panels
treePanel = uipanel(fig,'Title','Struct Tree','Position',[0.01 0.02 0.35 0.96]);
propPanel = uipanel(fig,'Title','Properties','Position',[0.37 0.02 0.62 0.96]);

% Buttons
uicontrol(treePanel,'Style','pushbutton','String','Refresh Tree',...
    'Units','normalized','Position',[0.02 0.01 0.45 0.06],...
    'Callback',@buildTree);

uicontrol(treePanel,'Style','pushbutton','String','Expand All',...
    'Units','normalized','Position',[0.52 0.01 0.46 0.06],...
    'Callback',@expandAll);

% Create tree UI
t = uitree(treePanel);
%t = uitree(treePanel,'Units','normalized', 'Position',[0.02 0.09 0.96 0.88]);
t.SelectionChangedFcn = @onNodeSelected;

% Initialize tree
buildTree();


% Wait for figure to close when function is called in blocking mode
uiwait(fig);

% Return updated struct
if isfield(data,'struct')
    updated = data.struct;
else
    updated = s;
end

% ----------------- Nested functions -----------------

    function buildTree(~,~)
        % Clear current tree
        delete(t.Children);
        root = uitreenode(t,'Text','root','NodeData',struct(), 'Icon',[]);
        % Recursively add nodes for each top-level field
        addStructNodes(root, data.struct, {});
        root.Expanded = true;
    end

    function addStructNodes(parentNode, st, path)
        % st is a struct; path is cell array of fieldnames leading here
        fn = fieldnames(st);
        for k=1:numel(fn)
            name = fn{k};
            value = st.(name);
            newPath = [path, {name}];
            if isstruct(value)
                % Node representing a struct: NodeData carries the path
                node = uitreenode(parentNode,'Text',name,'NodeData',newPath);
                node.Expanded = false;
                addStructNodes(node, value, newPath);
            else
                % Leaf node: show name and summary string
                summary = summarizeValue(value);
                txt = sprintf('%s: %s', name, summary);
                % store path and value in NodeData
                node = uitreenode(parentNode,'Text',txt,'NodeData',newPath);
            end
        end
    end

    function s = summarizeValue(v)
        % Create a short description for display in leaf text
        if isnumeric(v)
            if isscalar(v)
                s = num2str(v);
            else
                s = sprintf('numeric[%d]',numel(v));
            end
        elseif islogical(v)
            if isscalar(v)
                s = mat2str(v);
            else
                s = sprintf('logical[%d]',numel(v));
            end
        elseif ischar(v)
            s = ['''' v ''''];
        elseif isstring(v)
            s = ['"' char(v) '"'];
        elseif iscellstr(v)
            s = sprintf('cellstr[%d]', numel(v));
        elseif iscell(v)
            s = sprintf('cell[%d]', numel(v));
        else
            s = class(v);
        end
    end

    function onNodeSelected(src,evt)
        % When a node is selected, determine which struct-level it refers to
        % and populate property panel with editable controls for non-struct fields
        selected = evt.SelectedNodes;
        if isempty(selected)
            return;
        end
        node = selected(1);
        path = node.NodeData;
        % NodeData is either {} for root, or cell array of field names
        if isempty(path)
            % root selected -> show top-level fields
            st = data.struct;
            currentPath = {};
        else
            % Navigate to the struct that contains fields
            % We want to show fields at the level of path. If the node is
            % itself a struct node, show its fields. If it is a leaf field,
            % then show the fields of its parent struct.
            % Check whether this node corresponds to a struct in data.struct:
            try
                st = getStructAtPath(data.struct, path);
                currentPath = path;
            catch
                % If path points to a leaf value, show parent struct
                if numel(path)>1
                    parentPath = path(1:end-1);
                    st = getStructAtPath(data.struct, parentPath);
                    currentPath = parentPath;
                else
                    % top-level leaf: parent is root struct
                    st = data.struct;
                    currentPath = {};
                end
            end
        end
        populateProperties(st, currentPath);
    end

    function populateProperties(st, path)
        % Remove all children of propPanel then create controls
        delete(propPanel.Children);
        % Put a label for the current path
        if isempty(path)
            pathStr = 'root';
        else
            pathStr = strjoin(path,'.');
        end
        uicontrol(propPanel,'Style','text','String',['Path: ' pathStr],...
            'Units','normalized','Position',[0.02 0.95 0.96 0.04],'HorizontalAlignment','left');
        % Get fields that are NOT struct
        fn = fieldnames(st);
        % Filter out struct fields (we only edit non-struct fields here)
        editFields = {};
        for k=1:numel(fn)
            if ~isstruct(st.(fn{k}))
                editFields{end+1} = fn{k}; %#ok<AGROW>
            end
        end

        if isempty(editFields)
            uicontrol(propPanel,'Style','text','String','No editable fields at this level (struct contains only nested structs).',...
                'Units','normalized','Position',[0.02 0.5 0.96 0.4],'HorizontalAlignment','left');
            return;
        end

        % Layout variables
        n = numel(editFields);
        margin = 0.02;
        h = (0.9 - margin) / max(n,1); % vertical slot height
        startY = 0.02;

        % For each editable field, create label and appropriate control
        for i = 1:n
            fname = editFields{i};
            fval = st.(fname);
            ypos = startY + (n-i)*h;

            % Label
            uicontrol(propPanel,'Style','text','String',fname,...
                'Units','normalized','Position',[0.02 ypos+0.01 0.35 h-0.02],...
                'HorizontalAlignment','left');

            % Choose appropriate editor based on type
            if ischar(fval) || (isstring(fval) && isscalar(fval))
                % edit box for char / simple string
                ctrl = uicontrol(propPanel,'Style','edit','String',char(fval),...
                    'Units','normalized','Position',[0.39 ypos+0.01 0.59 h-0.02],...
                    'Callback',@(src,evt) onValueChanged(path, fname, get(src,'String')));
            elseif isscalar(fval) && isnumeric(fval)
                % numeric scalar
                ctrl = uicontrol(propPanel,'Style','edit','String',num2str(fval),...
                    'Units','normalized','Position',[0.39 ypos+0.01 0.59 h-0.02],...
                    'Callback',@(src,evt) onNumericChanged(path, fname, get(src,'String')));
            elseif islogical(fval) && isscalar(fval)
                % checkbox for logical scalar
                ctrl = uicontrol(propPanel,'Style','checkbox','Value',logical(fval),...
                    'Units','normalized','Position',[0.39 ypos+0.01 0.2 h-0.02],...
                    'Callback',@(src,evt) onBoolChanged(path, fname, get(src,'Value')));
            else
                % For vectors/cells/others show an edit button that opens a dialog
                uicontrol(propPanel,'Style','pushbutton','String','Edit Value',...
                    'Units','normalized','Position',[0.39 ypos+0.01 0.59 h-0.02],...
                    'Callback',@(~,~) onEditComplex(path, fname));
            end
        end

        % Small instruction
        uicontrol(propPanel,'Style','text','String','Edit values and they will be updated immediately in the struct. Close window to finish.',...
            'Units','normalized','Position',[0.02 0.91 0.96 0.04],'HorizontalAlignment','left','FontSize',9);
    end

    function onValueChanged(path, fname, newStr)
        % Generic change for string-like fields
        writeValueAtPath(path, fname, newStr);
        rebuildNodeText(path, fname);
    end

    function onNumericChanged(path, fname, str)
        % Try to parse numeric; if fail, show warning
        val = str2num(str); %#ok<ST2NM>
        if isempty(val)
            warndlg(sprintf('Cannot parse numeric from "%s". Value not changed.', str),'Invalid numeric','modal');
            return;
        end
        % If parsed matrix/vector, store as numeric array
        writeValueAtPath(path, fname, val);
        rebuildNodeText(path, fname);
    end

    function onBoolChanged(path, fname, val)
        writeValueAtPath(path, fname, logical(val));
        rebuildNodeText(path, fname);
    end

    function onEditComplex(path, fname)
        % Show a dialog to edit complex data (vectors, cells, etc.)
        cur = getValueAtPath(path, fname);
        answer = inputdlg({'Edit value (MATLAB expression):'},...
            sprintf('Edit %s', fname), [1 60], {mat2str(cur)});
        if isempty(answer)
            return;
        end
        try
            newVal = evalin('base', answer{1});
        catch
            try
                newVal = eval(answer{1});
            catch ME
                errordlg(['Error evaluating expression: ' ME.message],'Eval error','modal');
                return;
            end
        end
        writeValueAtPath(path, fname, newVal);
        rebuildNodeText(path, fname);
        % Refresh properties for the current selection
        % Find selected node and re-trigger selection
        sel = t.SelectedNodes;
        if ~isempty(sel)
            onNodeSelected(t, struct('SelectedNodes',sel));
        end
    end

    function v = getValueAtPath(path, fname)
        % Return the value of field fname at struct location path
        if isempty(path)
            parent = data.struct;
        else
            parent = getStructAtPath(data.struct, path);
        end
        v = parent.(fname);
    end

    function writeValueAtPath(path, fname, value)
        % Write back value into data.struct at location (path).(fname)
        if isempty(path)
            data.struct.(fname) = value;
        else
            % Use dynamic assignment by constructing a set command
            % Example: data.struct.a.b.c = value;
            % We implement with recursion to avoid eval when possible
            data.struct = setFieldAtPath(data.struct, path, fname, value);
        end
    end

    function S = setFieldAtPath(S, path, fname, value)
        % Set S.(path{1}).(path{2})... .(fname) = value
        if isempty(path)
            S.(fname) = value;
            return;
        end
        key = path{1};
        if numel(path)==1
            % modify this level
            sub = S.(key);
            sub.(fname) = value;
            S.(key) = sub;
        else
            sub = S.(key);
            sub = setFieldAtPath(sub, path(2:end), fname, value);
            S.(key) = sub;
        end
    end

    function st = getStructAtPath(S, path)
        % Return S.(path{1}).(path{2})... as struct. Throws error if not struct.
        st = S;
        for i=1:numel(path)
            fld = path{i};
            if ~isfield(st, fld)
                error('Field not found');
            end
            st = st.(fld);
            if ~isstruct(st) && i < numel(path)
                error('Not a struct at intermediate path');
            end
        end
        if ~isstruct(st)
            error('Target is not a struct');
        end
    end

    function rebuildNodeText(path, fname)
        % Update the tree node text for the leaf corresponding to path+fname
        % We traverse tree children to find matching node by comparing NodeData
        nodes = findAllNodes(t);
        targetPath = [path, {fname}];
        for ii=1:numel(nodes)
            nd = nodes(ii);
            if isequal(nd.NodeData, targetPath)
                % Update text: name: summary
                try
                    v = getValueAtPath(path, fname);
                catch
                    v = [];
                end
                nd.Text = sprintf('%s: %s', fname, summarizeValue(v));
                break;
            end
        end
    end

    function nodes = findAllNodes(treeRoot)
        % Return array of all nodes in the tree (breadth-first)
        nodes = [];
        stack = treeRoot.Children;
        while ~isempty(stack)
            nd = stack(1);
            stack(1) = [];
            nodes(end+1) = nd; %#ok<AGROW>
            if ~isempty(nd.Children)
                stack = [stack; nd.Children];
            end
        end
    end

    function expandAll(~,~)
        nodes = findAllNodes(t);
        for ii=1:numel(nodes)
            nodes(ii).Expanded = true;
        end
    end

    function onClose(~,~)
        % On close, resume and store struct
        uiresume(fig);
        delete(fig);
    end

end
