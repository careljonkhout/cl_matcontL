classdef EMat < handle
    %EMAT Embedded Matlab templating
    %
    %   EMat class provides a templating system in Matlab like Ruby's ERB
    %   system. Matlab code can be embedded inside any text document to
    %   easily control the document generation flow.
    %
    %   A simple example is illustrated here:
    %     >> x = 42;
    %     >> tmpl = '    The value of x is <%= x %>';
    %     >> obj = EMat(tmpl);
    %     >> disp(obj.render);
    %         The value of x is 42
    % 
    % Synopsis:
    %
    %   obj = EMat
    %   obj = EMat( S )
    %   obj = EMat( file_path )
    %
    %   obj.set( S )
    %   obj.set( file_path )
    %
    %   S = obj.render()
    %   obj.render( file_path )
    %
    %   obj = EMat creates an empty EMat object. EMat object accepts
    %   template string either by string variable S or by specifying a
    %   path to the template text file_path. To set a template to the
    %   object, use obj.set(S) or obj.set(file_path). EMat(S) and
    %   EMat(file_path) are the shorthand for obj=EMat; obj.set(...);
    %
    %   S = obj.render() returns a string of the rendered document.
    %   obj.render(file_path) instead renders output to a file specified
    %   by the file_path.
    %
    % Properties:
    %
    %     tmpl_path:  template file path. Use set() method to change
    %          tmpl:  template string. Use set() method to change
    %        errchk:  logical flag to enable/disable syntax check
    %                 (default: true)
    %          trim:  logical flag to enable/disable whitespace trim when
    %                 suppresseing newline at the end (default: true)
    %
    % Template format:
    %   
    %   Any text document can embed matlab code with the following syntax.
    %   
    %   <%  stmt  %> matlab statement
    %   <%  stmt -%> matlab statement with newline suppression at the end
    %   <%= expr  %> matlab expression with rendering
    %   <%# comt  %> comment line
    %   <%# comt -%> comment line with newline suppression at the end
    %   <%% %%>      escape to render '<%' or '%>', respectively
    %
    %   <%= expr %> renders output of the matlab expression to the output.
    %   Note that numeric variables will be converted to string by
    %   NUM2STR(). When -%> is specified at the end of the line in
    %   statement or comment, a following newline will be omitted from the
    %   rendering. Any other texts appearing outside of these special
    %   brackets are rendered as is. When trim property is set true,
    %   leading whitespace in the template is also removed from the output
    %   with newline suppression syntax.
    %
    % Example:
    %
    %   <!-- template.html.emat -->
    %   <html>
    %   <head><title><%= t %></title></head>
    %   <body>
    %   <%# this is a comment line -%>
    %   <p><%= a %></p>
    %   <ul>
    %   <% for i = 1:3 -%>
    %     <li><%= i %></li>
    %   <% end -%>
    %   </ul>
    %   </body>
    %   </html>
    %
    %   % In your matlab code
    %   % Prepare variables used in the template
    %   t = 'My template document';
    %   a = 10;
    %   
    %   % Create an EMat object
    %   obj = EMat('/path/to/template.html.emat');
    %   
    %   % Render to a file
    %   obj.render('/path/to/rendered.html');
    %
    % See also NUM2STR, FPRINTF

    % Revision 0.1   July 28, 2011
    % Revision 0.2   August 10, 2011
    % Revision 0.3   August 11, 2011
    % Revision 0.4   August 12, 2011
    % Revision 0.5   December 22, 2011
    %
    % You can redistribute this software under BSD license
    % Copyright (c) 2011 Kota Yamaguchi
    
    properties (Constant, Hidden)
        SPLIT_REGEXP = '<%%|%%>|<%=|<%#|<%|-%>|%>|\n'
    end
    
    properties (SetAccess = protected, Hidden)
        stag = ''
        last = ''
        content = ''
        stmts = {}
        script = ''
        out = ''
        fid = 1       % Default: 1=stdout
    end
    
    properties (SetAccess = protected)
        tmpl_path = ''
        tmpl = ''
    end
    
    properties
        errchk = true % Default: true
        trim   = true % Default: true
    end
    
    methods
        function [ obj ] = EMat(input)
            %EMat create a new EMat object
            if nargin > 0, obj.set(input); end
        end
        
        function [] = set(obj, input)
            %SET set template to render
            
            % Check input args
            if ischar(input)
                if exist(input,'file')
                    % Load input file if a pathname specified
                    obj.tmpl_path = input;
                    f = fopen(input,'r');
                    obj.tmpl = fscanf(f,'%c',inf);
                    fclose(f);
                else
                    obj.tmpl_path = '';
                    obj.tmpl = input;
                end
            else
                error('EMat:set:invalidInput',...
                    'Input argument must be a path or a string');
            end
        end
        
        function [ s ] = render(obj, output)
            %RENDER render template
            
            % Set fid if optional output path specified
            if nargin > 1 && nargout > 0
                warning('EMat:render:unsupported',...
                    ['output to string is unsupported '...
                    'when exporting to a file']);
            end
            if nargin > 1 && ischar(output)
                obj.fid = fopen(output,'w');
            end
            
            try
                % Compile template
                src = obj.compile(obj.tmpl);
                % Check syntax error
                if obj.errchk, obj.syntax_check(src); end
                % Render template in the caller context
                s = evalc('evalin(''caller'',src);');
            catch e
                if nargin > 1 && ischar(output)
                    fclose(obj.fid);
                    obj.fid = 1;
                end
                rethrow(e);
            end
            
            if nargin > 1 && ischar(output)
                fclose(obj.fid);
                obj.fid = 1;
            end
        end
        
        function [] = set.errchk(obj, value)
            %SET.ERRCHK set errchk flag
            obj.errchk = logical(value);
        end
        
        function [] = set.trim(obj, value)
            %SET.TRIM set trim flag
            obj.trim = logical(value);
        end
    end
    
    methods (Access = protected)
        function [ s ] = compile(obj, s)
            %COMPILE compile template string
            obj.scan(s);
            if ~isempty(obj.content), obj.push_print; end
            if ~isempty(obj.stmts), obj.cr; end
            s = obj.script;
            obj.clean;
        end
        
        function [] = scan(obj, s)
            %SCAN scan and tokenize text
            [match,lines] = regexp(s,'\n','match','split');
            lines = strcat(lines, [match,{''}]);
            for i = 1:length(lines)
                line = lines{i};
                [match,tokens] = regexp(line,EMat.SPLIT_REGEXP,...
                    'match','split');
                tokens = [(tokens); [match,{''}]];
                tokens = tokens(:);
                tokens(cellfun(@isempty,tokens)) = [];
                if obj.trim
                    % Trim whitespace if the end is '-%>\n'
                    if numel(tokens)>2 &&...
                       strcmp(tokens{end},char(10)) && ...
                       strcmp(tokens{end-1},'-%>')
                        ind = 1:find(strcmp(tokens,'<%'),1)-1;
                        tokens(ind) = strtrim(tokens(ind));
                    end
                    tokens(cellfun(@isempty,tokens)) = [];
                end
                for j = 1:numel(tokens)
                    obj.process(tokens{j});
                end
            end
        end
        
        function [] = process(obj, tok)
            %PROCESS parse tokens
            if isempty(obj.stag) % State 1: stag doesn't exist
                switch tok
                    case {'<%', '<%=', '<%#'}
                        obj.stag = tok;
                        if ~isempty(obj.content), obj.push_print; end
                        obj.content = '';
                    case 10 % '\n'
                        if ~strcmp(obj.last,'-%>')
                            obj.content = [obj.content,10];
                        end
                        obj.push_print;
                        obj.cr;
                        obj.content = '';
                    case '<%%'
                        obj.content = [obj.content,'<%%'];
                    otherwise
                        obj.content = [obj.content,tok];
                end
            else % State 2: stag exists
                switch tok
                    case {'%>','-%>'}
                        switch obj.stag
                            case '<%'
                                if obj.content(end)==10 % '\n'
                                    obj.content(end) = [];
                                    obj.push;
                                    obj.cr;
                                else
                                    obj.push;
                                end
                            case '<%='
                                obj.push_insert;
                            case '<%#'
                                % do nothing
                        end
                        obj.stag = '';
                        obj.content = '';
                    case '%%>'
                        obj.content = [obj.content,'%%>'];
                    otherwise
                        obj.content = [obj.content,tok];
                end
            end
            obj.last = tok;
        end
        
        function [] = push(obj)
            %PUSH add raw stmt
            obj.stmts = [obj.stmts,{obj.content}];
        end
        
        function [] = push_print(obj)
            %PUSH_PRINT add print stmt
            if ~isempty(obj.content)
                obj.stmts = [obj.stmts,...
                    {['fprintf(',num2str(obj.fid),',',...
                    obj.dump(obj.content),')']}];
            end
        end
        
        function [] = push_insert(obj)
            %PUSH_INSERT add insertion stmt
            obj.stmts = [obj.stmts,...
                {['fprintf(',num2str(obj.fid),...
                ',EMat.str(',obj.content,'))']}];
        end
        
        function [] = cr(obj)
            %CR flush stmts to script
            s = strtrim(obj.stmts);
            delim = repmat({';'},1,length(obj.stmts));
            % mlint complains about semi-colon after else stmt
            for i = find(strcmp(s,'else')), delim{i} = ' '; end
            s = [s;delim];
            obj.script = [obj.script, [s{:}]];
            obj.stmts = {};
            obj.script = sprintf('%s\n',obj.script);
        end
        
        function [] = clean(obj)
            %CLEAN reset properties
            obj.stag = '';
            obj.content = '';
            obj.stmts = {};
            obj.script = '';
        end
        
        function [ inform ] = syntax_check(obj, src)
            %SYNTAX_CHECK export src into tempfile and check error
            
            % Escape chars
            src = regexprep(src,'\\','\\\\');
            src = regexprep(src,'%','%%');
            
            % Write to tempfile
            file_path = [tempname,'.m'];
            f = fopen(file_path,'w');
            fprintf(f,src);
            fclose(f);
            
            % Check syntax
            inform = mlint(file_path,'-struct');
            if ~isempty(inform)
                warning('EMat:syntax_check:syntaxWarning','');
                for i = 1:length(inform)
                    fprintf('%s:line %d: %s\n',...
                        obj.tmpl_path, inform(i).line, inform(i).message);
                end
            end
            
            % Delete tempfile
            delete(file_path);
        end
    end
    
    methods(Static)        
        function [ s ] = dump(s)
            %DUMP escape and enclose char
            s = regexprep(s,'''','''''');
            s = regexprep(s,'\\','\\\\');
            s = strrep(s,'%','%%');
            s = strrep(s,char(10),'\n');
            s = ['''',s,''''];
        end
        
        function [ s ] = str(s)
            %STR string conversion
            if isnumeric(s)
                s = num2str(s);
            else
                s = char(s);
            end
        end
    end
    
end

