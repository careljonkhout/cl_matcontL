classdef System_of_ODEs < matlab.mixin.CustomDisplay
  % 
  % see README for instructions on how to use this tool to generate
  % system files for (cl_)matcont/matcontL
  %
  % defining properties
  % These properties are supplied by the user and fully define the system.
  properties (SetAccess=immutable)
    name              % name of the system, used to create a filename
                      % char array
    input_vars_str    % list of variables of the system
                      % coded as a char array separated by spaces
    input_pars_str    % list of parameters of the system
                      % coded as a char array separated by spaces
    time              % name of the variable that represents time
                      % char array
    max_ord_derivatives
                      % the maximum order of derivatives that is 
                      % to be computed
                      % integer
    rhs               % an n by 1 cell array of character vectors
    
    output_type       % string. should equal odefile, C, or odemex
    
    app               % optional: GUI class that calls System_of_ODEs
                      % if app is given then status notifications
                      % will be passed to app by
                      % calling app.updatStatus(status)
                      
  end
  
  % These properties are computed from the defining propeties when an
  % instance of this class is created.
  properties (SetAccess=immutable)
    input_vars % cell array of char arrays
    input_pars % cell array of char arrays
    
    % symbols (i.e. variables and parameters) supllied by the user 
    % cell array of char arrays
    % the order corresponds to internal_syms and output_syms
    % i.e. input_syms{i} corresponds to internal_syms{i}
    % and output_syms{i}
    input_syms                 
    
    % symbols used by this class when computing derivatives.
    % cell array of char arrays
    internal_syms
    internal_vars % cell array of char arrays
    internal_pars % cell array of char arrays
    
    internal_vars_str  % char array
    internal_pars_str % char array
    
    % symbols in the form that is written to the system file
    % cell array of char arrays 
    output_syms
    
    % all fields below contain char arrays
    syms_arg   
    formatted_rhs
    parameter_arguments
    jacobian                          = '[]'
    jacobian_handle                   = '[]'
    jacobian_params                   = '[]'
    jacobian_params_handle            = '[]'
    hessians                          = '[]'
    hessians_handle                   = '[]'
    hessians_params                   = '[]'
    hessians_params_handle            = '[]'  
    third_order_derivatives           = '[]'
    third_order_derivatives_handle    = '[]' 
    fourth_order_derivatives          = '[]'
    fourth_order_derivatives_handle   = '[]'
    fifth_order_derivatives           = '[]'
    fifth_order_derivatives_handle    = '[]'
    
  end

  
  methods (Access = protected)
    % this method defines the customized output that is displayed when
    % inspecting a System_of_ODEs on the command line or when a System_of_ODEs
    % is printed using the disp function see:
    % https://mathworks.com/help/matlab/ref/matlab.mixin.util.propertygroup-class.html
    function propgrp = getPropertyGroups(s)
      my_rhs = replace_symbols(strjoin(s.rhs, ', '), ...
        s.internal_syms, s.input_syms);
      props = struct( ...
         'name',                          s.name, ...
         'variables',                     s.input_vars_str, ...
         'parameters',                    s.input_pars_str, ...
         'maximum_order_of_derivatives',  num2str(s.max_ord_derivatives), ...
         'time_variable',                 s.time, ...
         'right_hand_side',               my_rhs ...
         );
      if s.max_ord_derivatives >= 1
        props.jacobian = replace_symbols(s.jacobian, ...
          s.output_syms, s.input_syms);
        props.jacobian_params = ...
          replace_symbols(s.jacobian_params, ...
            s.output_syms, s.input_syms);
      end
      if s.max_ord_derivatives >= 2
        props.hessians = replace_symbols(s.hessians, ...
          s.output_syms, s.input_syms);
        props.hessians_params = ...
          replace_symbols(s.hessians_params, ...
            s.output_syms, s.input_syms);
      end
      if s.max_ord_derivatives >= 3
        props.third_order_derivatives = replace_symbols(...
          s.third_order_derivatives, s.output_syms, s.input_syms);
      end
      if s.max_ord_derivatives >= 4
        props.fourth_order_derivatives = replace_symbols(...
          s.fourth_order_derivatives, ...
          s.output_syms, s.input_syms);
      end
      if s.max_ord_derivatives >= 5
        props.fifth_order_derivatives = replace_symbols(...
          s.fifth_order_derivatives, ...
          s.output_syms, s.input_syms);
      end

      propgrp = matlab.mixin.util.PropertyGroup(props);
    end
  end

  methods
    % constructor
    % rhs may be supplied as either a n by 1 string
    % a 1 by n string array
    % or a n by m character array, where 
    % n is the number of equations
    % and m is the length of the longest equation
    function s = System_of_ODEs(...
        name,...
        variables_str,...
        parameters_str,...
        time,...
        max_ord_derivatives,...
        rhs,...
        app,...
        output_type)
            
      % set defining properties
      s.name                = strtrim(char(name));
      s.input_vars_str      = strtrim(char(variables_str));
      s.input_pars_str      = strtrim(char(parameters_str));
      s.rhs                 = strtrim(cellstr(rhs));
      s.max_ord_derivatives = max_ord_derivatives;
      s.time                = strtrim(char(time));
      s.app                 = app;
      s.output_type         = output_type;
      
      
            
      % set derived properties
      s.input_vars = regexp(s.input_vars_str, '( |,)+','split');    
      s.input_pars = regexp(s.input_pars_str, '( |,)+','split');
      if isempty(s.input_pars{1})
        s.input_pars = {};
        s.input_pars_str = '';
        s.parameter_arguments = '';
      else
        s.parameter_arguments = s.generate_parameter_arguments();
        s.input_pars_str = strjoin(s.input_pars, ' ');
      end
      s.input_vars_str  = strjoin(s.input_vars, ' ');
      [s.input_syms, s.internal_syms, s.output_syms] = s.generate_sym_lists;
      
      s.verify_inputs()

      for i=1:length(s.rhs)
        s.rhs{i} = replace_symbols(s.rhs{i}, s.input_syms, s.internal_syms);
      end
  
      s.internal_vars_str = s.input_to_internal(s.input_vars_str);
      s.internal_vars = strsplit(s.internal_vars_str);
      
      if ~isempty(s.input_pars)
        s.internal_pars_str = s.input_to_internal(s.input_pars_str);
        s.internal_pars = strsplit(s.internal_pars_str);
      end      
      
      s.syms_arg = [s.internal_vars_str ' ' s.internal_pars_str];
        
      s.formatted_rhs = s.format_rhs();
      
      if (s.max_ord_derivatives >= 1)
        s.jacobian                = s.compute_jacobian_for_variables();
        s.jacobian_params = s.compute_jacobian_params();
        switch s.output_type
          case 'c'
            s.jacobian_handle         = sprintf('@%s_jacobian'       , s.name);
            s.jacobian_params_handle  = sprintf('@%s_jacobian_params', s.name);
          case 'odefile'
            s.jacobian_handle         = '@jacobian';
            s.jacobian_params_handle  = '@jacobian_params';
          case 'odemex'
        end
      end
      if (s.max_ord_derivatives >= 2)
        s.show_status('Computing second order derivatives...')
        s.hessians                = s.compute_hessians('vars');
        s.hessians_params         = s.compute_hessians('pars');
        s.hessians_handle         = '@hessians';
        s.hessians_params_handle  = '@hessians_params';
      end
      if (s.max_ord_derivatives >= 3)
        s.show_status('Computing third order derivatives...')
        s.third_order_derivatives = s.compute_3rd_ord_derivatives;
        s.third_order_derivatives_handle = '@third_order_derivatives';
      end
      if (s.max_ord_derivatives >= 4)
        s.show_status('Computing fourth order derivatives...')
        s.fourth_order_derivatives = s.compute_4th_ord_derivatives;
        s.fourth_order_derivatives_handle = '@fourth_order_derivatives';
      end
      if (s.max_ord_derivatives >= 5)
        s.show_status('Computing fifth order derivatives...')
        s.fifth_order_derivatives = s.compute_5th_ord_derivatives;
        s.fifth_order_derivatives_handle = '@fifth_order_derivatives';
      end
    end
    
    function out = input_to_internal(s, in)
      out = replace_symbols(in, s.input_syms, s.internal_syms);
    end
    
    function generate_file(s)
      switch s.output_type
        case 'C'
          s.generate_c_files();
        case 'odemex'
          s.generate_odemex_files();
        case 'odefile'
          path = get_path();

          filename = fullfile(path, '../Systems', [s.name,'.m']);
          file_contents = emat2('system.m.emat');
          fileID = fopen(filename,'w');
          fprintf(fileID,'%s',file_contents);
      end
    end
    
    function generate_c_files(s)
      path = get_path();
      
      filename = fullfile(path, '../Systems', [s.name '_mex.m']);
      file_contents = emat2('system_mex.m.emat');
      fileID = fopen(filename,'w');
      fprintf(fileID,'%s',file_contents);
      
      filename = fullfile(path, '../Systems', [s.name '_dydt.c']);
      file_contents = emat2('system_dydt.c.emat');
      fileID = fopen(filename,'w');
      fprintf(fileID,'%s',file_contents);
      
      mex_filename = fullfile(path, '../Systems', [s.name '_dydt']);
      mex(filename,'-output', mex_filename);
      
      if (s.max_ord_derivatives >= 1)
        filename = fullfile(path, '../Systems', [s.name '_jacobian.c']);
        file_contents = emat2('system_jacobian.c.emat');
        fileID = fopen(filename,'w');
        fprintf(fileID,'%s',file_contents);
        
        mex_filename = fullfile(path, '../Systems', [s.name '_jacobian']);
        mex(filename,'-output', mex_filename);
        
        filename = fullfile(path, '../Systems', [s.name '_jacobian_params.c']);
        file_contents = emat2('system_jacobian_params.c.emat');
        fileID = fopen(filename,'w');
        fprintf(fileID,'%s',file_contents);
        
        mex_filename = fullfile(path, '../Systems',[s.name '_jacobian_params']);
        mex(filename,'-output', mex_filename);
      end
    end
    
    function generate_odemex_files(s)
      path = get_path();
      filename = fullfile(path, '../Systems', [s.name '_odemex_dydt.m']);
      file_contents = emat2('odemex_dydt.m.emat');
      fileID = fopen(filename,'w');
      fprintf(fileID, '%s', file_contents);
    end
    
    function verify_inputs(s)
      all_equations = join(s.rhs);
      if isempty(s.name)
        s.throwException('You must provide a system name.');
      elseif isempty(s.input_vars_str)
        s.throwException('You must provide at least one variable.');
      elseif isempty(s.rhs)
        s.throwException('You must specify a right hand side.');
      elseif length(s.rhs) ~= length(s.input_vars)
        msg = 'The number of variables must equal the number of equations.';
        s.throwException(msg);
      elseif contains(all_equations, ';')
        s.throwException('The equations may not contain a semicolon.');
      elseif contains(all_equations, '=')
        s.throwException('The equations may not contain an equals sign.');
      elseif strcmp(s.output_type, 'C') && (s.max_ord_derivatives > 1)
        s.throwException( ...
         'C-output of derivatives of order greater than 1 is not implemented.');
      elseif ~ any(strcmp({'odefile', 'C', 'odemex'},s.output_type))
        s.throwException([s.output_type ' is not a valid output type']);
      else
        for i=1:length(s.rhs)
          if ~ System_of_ODEs.match_parentheses(s.rhs{i})
            s.throwException(['There is a mismatched parenthesis in' ...
             ' right hand side expression number ' int2str(i) '.']);
          end
        end
      end
    end

    function str = generate_parameter_arguments(s)
      new_parameters = strcat('par_', s.input_pars);
      str = strjoin(new_parameters, ', ');
    end

    function formatted_rhs = format_rhs(s)
      switch s.output_type
        case 'C'
          eval(['syms ' s.syms_arg]);
          dydt_elements = cell(length(s.rhs),1);
          for i = 1 : length(s.rhs)

            raw_c_code = eval(sprintf('ccode(%s)',s.rhs{i}));
            tokens    = regexp(raw_c_code, '.*=(.*);', 'tokens');
            c_code    = tokens{1}{1};
            dydt_elements{i} = sprintf('  dydt[%d] = %s;', i-1, c_code);
          end
          formatted_rhs = strjoin(dydt_elements, '\n');
        case 'odemex'
          eval(['syms ' s.syms_arg]);
          dydt_elements = cell(length(s.rhs),1);
          for i = 1 : length(s.rhs)
            raw_odemex_code = eval(sprintf('ccode(%s)',s.rhs{i}));
            tokens    = regexp(raw_odemex_code, '.*=(.*);', 'tokens');
            odemex_code    = tokens{1}{1};
            dydt_elements{i} = sprintf('dx(%d) = %s;', i, odemex_code);
          end
          formatted_rhs = strjoin(dydt_elements, '\n');
        case 'odefile'
          formatted_rhs = ['[' strjoin(s.rhs, '; ') ']'];
      end
      formatted_rhs = replace_symbols(formatted_rhs, ...
                                                s.internal_syms, s.output_syms);
    end

    function [in, intermediate, out] = generate_sym_lists(s)
      vars_len = length(s.input_vars);
      symbols_length = vars_len + length(s.input_pars);
      in           = cell(1, symbols_length);
      intermediate = cell(1, symbols_length);
      out          = cell(1, symbols_length);
      for i=1:length(s.input_vars)
        in{i}           = s.input_vars{i};
        intermediate{i} = ['v_', s.input_vars{i}];
        switch s.output_type
          case 'C'
          out{i}          = sprintf('y[%d]', i-1);
          case 'odemex'
          out{i}          = sprintf('x(%d)', i  );
          case 'odefile'
          out{i}          = sprintf('y(%d)', i  );
        end
      end
      for i = 1:length(s.input_pars)
        in{vars_len + i}           = s.input_pars{i};
        intermediate{vars_len + i} = ['p_', s.input_pars{i}];
        switch s.output_type
          case 'C'
          out{vars_len + i}          = sprintf('parameters[%d]', i - 1);
          case 'odemex'
          out{vars_len + i}          = sprintf('p(%d)'         , i    );
          case 'odefile'
          out{vars_len + i}          = ['par_', s.input_pars{i}]; 
        end
      end
    end

    function jacobian = compute_jacobian_safe(s, input_vars_str)
      eval_str = sprintf('jacobian([%s],[%s])', ...
          strjoin(s.rhs, ','), input_vars_str);
      jacobian = s.safe_eval_w_syms(eval_str);
      switch s.output_type
        
        case 'C'
          
        c_code_of_jac_elements = s.to_c_code(jacobian);
        jacobian_elements      = cell(numel(jacobian),1);
        j = 1;
        for i = 1 : numel(jacobian)
          if ~ isempty(c_code_of_jac_elements{i})
            jacobian_elements{j} = sprintf('  jac[%d] = %s;', i-1, ...
                                                     c_code_of_jac_elements{i});
	          j = j + 1;
          end
        end
        number_of_nonzeros = j;
        % todo: remove empty strings
        jacobian = strjoin(jacobian_elements{1:number_of_nonzeros}, '\n');
        
        case 'odemex'
          
        jacobian = s.symbolic_mat_to_matlab_code(jacobian);
        
        case 'odefile'
          
        jacobian = s.symbolic_mat_to_matlab_code(jacobian);
       
      end
      jacobian = replace_symbols(jacobian, s.internal_syms, s.output_syms);
    end
    
    function symbolic_mat = symbolic_mat_to_matlab_code(s, symbolic_mat)
      eval(['syms ' s.syms_arg]);
      symbolic_mat = char(matlabFunction(symbolic_mat));
      % jacobian now is a char array with the matlab code of the Jacobian
      % matrix. For instance if the jacobian is [p_a ; p_b], then the variable
      % "jacobian" now contains the string "@(p_a, p_b) [p_a, p_b]". If the
      % number of dimensions of the Jacobian is greater than 1, then the
      % variable "jacobian" contains a substring 'reshape'.
      if contains(symbolic_mat, 'reshape')
        tokens   = regexp(symbolic_mat, '(reshape.*)', 'tokens');
        symbolic_mat = tokens{1}{1};
      else
        tokens   = regexp(symbolic_mat, '(\[.*\])', 'tokens');
        symbolic_mat = tokens{1}{1};
      end
    end
    
    function c_code = to_c_code(s, expressions)
      eval(['syms ' s.syms_arg]);
      c_code = cell(numel(expressions),1);
      for i = 1 : numel(expressions)
        if ~ strcmp(char(expressions(i)),'0')
          my_c_code = ccode(expressions(i));
          tokens    = regexp(my_c_code, '.*=(.*);', 'tokens');
          c_code{i} = tokens{1}{1};
        end
      end
    end
    
    function matlab_code = to_matlab_code(s, expressions)
      eval(['syms ' s.syms_arg]);
      matlab_code = cell(numel(expressions),1);
      for i = 1 : numel(expressions)
        %if ~ strcmp(char(expressions(i)),'0')
          my_matlab_code = char(matlabFunction(expressions(i)));
          tokens    = regexp(my_matlab_code, '.*=(.*);', 'tokens');
          matlab_code{i} = tokens{1}{1};
        %end
      end
    end
        
    % returns the jacobian matrix as a string of matlab code
    % with the i-th variable replaced with y(i)
    % and a parameter a replaced with par_a
    function jacobian = compute_jacobian_for_variables(s)
      jacobian = s.compute_jacobian_safe(s.internal_vars_str);
    end

    function jacobian = compute_jacobian_params(s)
      jacobian = s.compute_jacobian_safe(s.internal_pars_str);
    end

    function d = compute_hessians(s, vars_or_pars)
      eval(['syms ' s.syms_arg]);
      switch vars_or_pars
        case 'vars'
        vars = s.internal_vars;
        case 'pars'
        vars = s.internal_pars;
      end
      len   = length(vars);
      d_sym = sym(zeros(length(s.rhs),len,len));
      for i = 1:length(s.rhs)
        for j = 1:length(vars)
          for k = 1:j
            evalstr = sprintf('diff(%s, %s, %s)', s.rhs{i}, vars{j}, vars{k});
            pd = eval(evalstr);
            d_sym(i,j,k) = pd;
            d_sym(i,k,j) = pd;
          end
        end
      end
      d = s.symbolic_mat_to_matlab_code(d_sym);
      d = replace_symbols(d, s.internal_syms, s.output_syms);
    end
    
    function d = compute_3rd_ord_derivatives(s)
      eval(['syms ' s.syms_arg]);
      variables = s.internal_vars;
      len   = length(variables);
      d_sym = sym(zeros(length(s.rhs),len,len,len));
      for i = 1:length(s.rhs)
        for j = 1:length(variables)
          for k = 1:j
            for l = 1:k
              evalstr = sprintf('diff(%s, %s, %s, %s)', ...
                  s.rhs{i}, variables{j}, ...
                  variables{k}, variables{l});
              pd = eval(evalstr);
              ps = perms([j,k,l]);
              for m =1:length(ps)
                d_sym(i,ps(m,1), ps(m,2), ps(m,3)) = pd;
              end
            end
          end
        end
      end
      d = s.symbolic_mat_to_matlab_code(d_sym);
      d = replace_symbols(d, s.internal_syms, s.output_syms);
    end

    function d = compute_4th_ord_derivatives(s)
      eval(['syms ' s.syms_arg]);
      variables = s.internal_vars;
      len = length(variables);
      d_sym = sym(zeros(length(s.rhs),len,len,len,len));
      for i = 1:length(s.rhs)
        for j = 1:length(variables)
          for k = 1:j
            for l = 1:k
              for m = 1:l
                evalstr = sprintf('diff(%s, %s, %s, %s, %s)', ...
                    s.rhs{i}, variables{j}, variables{k}, variables{l}, ...
                    variables{m});
                pd = eval(evalstr);
                p = perms([j, k, l, m]);
                for pi =1:length(p)
                  d_sym(i,p(pi, 1), p(pi, 2), p(pi, 3), p(pi, 4)) = pd;
                end
              end
            end
          end
        end
      end
      d = s.symbolic_mat_to_matlab_code(d_sym);
      d = replace_symbols(d, s.internal_syms, s.output_syms);
    end

    function d = compute_5th_ord_derivatives(s)
      eval(['syms ' s.syms_arg]);
      variables = s.internal_vars;
      len = length(variables);
      d_sym = sym(zeros(length(s.rhs),len,len,len,len,len));
      for i = 1:length(s.rhs)
        for j = 1:length(variables)
          evalstr = sprintf('diff(%s, %s)', s.rhs{i}, variables{j});
          % the pd1, pd2, pd3 will be used via the "eval" function
          pd1 = eval(evalstr); %#ok<NASGU>
          for k = 1:j
            evalstr = sprintf('diff(pd1, %s)', variables{k});
            pd2 = eval(evalstr); %#ok<NASGU>
            for l = 1:k
              evalstr = sprintf('diff(pd2, %s)', variables{l});
              pd3 = eval(evalstr); %#ok<NASGU>
              for m = 1:l
                evalstr = sprintf('diff(pd3, %s)', variables{m});
                pd4 = eval(evalstr);
                for n = 1:m
                  if (pd4 == 0)
                    pd = sym(0);
                  else
                    evalstr = sprintf('diff(pd4, %s)',variables{n});
                    pd = eval(evalstr);
                  end
                  p = perms([j, k, l, m, n]);
                  for pi =1:length(p)
                    d_sym(i,p(pi,1), p(pi,2), p(pi,3), p(pi,4), p(pi,5)) = pd;
                  end
                end
              end
            end
          end
        end
      end
      d = s.symbolic_mat_to_matlab_code(d_sym);
      d = replace_symbols(d, s.internal_syms, s.output_syms);
    end

    function result=safe_eval_w_syms(s,evalstr)
      try 
        eval(['syms ' s.syms_arg]);
        result = eval(evalstr);        
      catch exception
        switch exception.identifier
          case'MATLAB:UndefinedFunction'
            regex = '''\w*''';
            variable = regexp(exception.message, regex, 'match');
            msg = 'The system equations contain the variable';
            msg = [msg ' ' variable{1} ' '];
            msg = [msg 'which you did not specify in the variables '];
            msg = [msg 'or the parameters. '];
            msg = [msg 'Add ' variable{1} ' to the variables '];
            msg = [msg 'or parameters and try again.'];
            MException('System_of_ODEs:UndefinedFunctionInEval', ...
              msg).throw;
%         case 'symbolic:sym:SymsCannotOvershadow'
%           exception.rethrow();
          otherwise
            msg = 'An error occured while ' ;
            msg = [msg 'evaluating the statement: " '];
            msg = [msg char(evalstr) ' "\n'];
            msg = [msg 'The error was: \" '];
            msg = [msg exception.message ' "\n'];
            msg = [msg 'The symbolic variables in this context were: '];
            msg = [msg char(s.syms_arg)];
            MException('System_of_ODEs:EvalError', msg).throw;
        end
      end
    end
    
    function show_status(s,status)
      if ~isempty(s.app)
        s.app.updateStatus(status);
      else
        disp(status)
      end
    end
     
  end
    
  methods(Static)

    
    function parentheses_match = match_parentheses(str)
      parentheses_level = 0;
      for i=1:length(str)
        if strcmp(str(i),'(')
          parentheses_level = parentheses_level + 1;
        elseif strcmp(str(i),')')
          parentheses_level = parentheses_level - 1;
          if (parentheses_level < 0)
            parentheses_match = false;
            return;
          end
        end
      end
      parentheses_match = parentheses_level == 0;
    end
        
    function throwException(msg)
      exception = MException('System_of_ODEs:BadInput', msg);
      throw(exception);
    end
    
    function s=new(name, var_str, par_str, time, max_ord, rhs, output_type)
       if nargin == 6
        output_type = 'odefile';
      end
      s = System_of_ODEs(name,var_str,par_str,time,max_ord,rhs,[],output_type);
    end
    
  end   
  
end