classdef SystemFileGenerator < matlab.mixin.CustomDisplay & handle
  %
  % written by Carel Jonkhout 2018-2020
  %
  % see Readme for instructions on how to use this tool to generate
  % system files for cl_matcontL
  %
  % defining properties
  % These properties are supplied by the user and fully define the system.
  properties
    name              % char array
                      % name of the system, used for the filename
                      
    input_pars_str    % list of parameters of the system
                      % coded as a char array separated by spaces
    time              % name of the variable that represents time
                      % char array
    max_ord_derivatives
                      % the maximum order of derivatives that is 
                      % to be computed
                      % integer
	 
    equations         % a 1d cell array of character vectors
                      % containing at least 1 differential equation
                      % and optionally partial expressions
                      % for example:
                      % 
                      % f = x + y  (optional partial expression)
                      % x' = f * y (differential equation)
                      % y' = f * x (differential equation)
                          
    
    output_type       % A char array. Should equal 'odefile', 'C', 'cvode'
                      % 'C_sparse', or 'cvode_sparse'
    
    app               % optional: GUI class that calls SystemFileGenerator
                      % if app is specified then status notifications
                      % will be passed to app by
                      % calling app.updatStatus(status)
                      
  end
  
  % The properties below are computed from the defining propeties when an
  % instance of this class is created.
  % The example value of each property is indicated for the input:
  % input_parameters = 'p q';
  % equations = { 'f = p^2'
  %               'g = q^2'
  %               'x'' = f*x - y'
  %               'y'' = x + g*y' }
  % output = 'odefile'
  
  
  properties
    input_vars           % cell array of char arrays. e.g. {'x', 'y'}
    input_pars           % cell array of char arrays. e.g. {'p', 'q'}
    input_intermediates  % cell array of char arrays. e.g. {'f', 'g'}
    
    % symbols (i.e. variables and parameters, and intermediate variables)
    % suplied by the user. The propery input_syms is cell array of char arrays
    % the order corresponds to internal_syms and output_syms i.e. input_syms{i}
    % corresponds to internal_syms{i} and output_syms{i}
    input_syms           % e.g. { 'x', 'y', 'p', 'q', 'f', 'g' }
    input_vars_str       % char array. e.g. 'x y'
    
    % the right hand sides of the differential equations stored as a cell array
    % of char arrays, with variables written in internal format
    % The purpose of the internal variable notation is to prevent local
    % variables in functions being overwritten by the symbols. 
    rhs                  % e.g. { 'i_f * v_x - v_y', 'v_x + i_g * v_y' }
    
    % the intermediate expressions stored as a cell array
    % of char arrays, with variables written in internal format
    intermediate_expressions % e.g. { 'i_f = p_p^2', i_g = p_q^2 }
    
    % the right hand sides of the system of ODEs with the intermediate
    % expressions substituted:
    rhs_with_substitutions  
    % e.g. { 'p_p^2 * v_x - v_y', 'v_x + p_q^2 * v_y' }
    
    % symbols used by this class when computing derivatives.
    % cell array of char arrays
    internal_syms % e.g. { 'v_x', 'v_y', 'p_p', 'p_q', 'i_f', 'i_g' }
    internal_vars % cell array of char arrays. e.g. {'v_x', 'v_y'} 
    internal_pars % cell array of char arrays. e.g. {'p_p', 'p_q'}
    
    internal_vars_str % char array. e.g. 'v_x v_y'
    internal_pars_str % char array. e.g. 'p_p p_q'
    
    % symbols in the form that is written to the system file
    % cell array of char arrays 
    output_syms % cell array of char arrays. e.g. 
    % if output is 'odefile':
    %      { 'y(1)', 'y(2)', 'par_p', 'par_q', 'i_f', 'i_g' } 
    % if output is 'C' or 'cvode'
    %      { 'y[0]', 'y[1]', 'parameters[0]', 'parameters[1]', 'i_f', 'i_g' } 
    
    
    jacobian_vars_sym 
    % Jacobian matrix of the ODE rhs w.r.t. the depedent variables 
    % stored as Matlab symbolic toolbox symbolic expressions
    
    jacobian_pars_sym
    % Jacobian matrix of the ODE rhs w.r.t. the parmeters
    % stored as Matlab symbolic toolbox symbolic expressions
    
    
    % all fields below contain char arrays
    syms_arg                          % e.g. 'v_x v_y p_p p_q i_f i_g'
    formatted_rhs                     % e.g.
    % if output is 'odefile'
    %  'i_f       = par_p^2 
    %   i_g       = par_q^2'
    %   dydt = [
    %       i_f*x - y(1)
    %       x + i_g*y(2)]; 
    parameter_arguments % e.g. par_p, par_q
    jacobian                          = '[]' 
    % if max_ord is zero jacobian will equal '[]'
    % if max_ord is greater than or equal to 1
    %   if  output_type is equal to 'odefile'
    %     'reshape([par_p^2,1], [1,par_q^2],2,2)'
    %   if output_type is equal to 'C' or 'cvode'
    %     JAC(0,0) = parameters[0]^2;
    %     JAC(1,0) = 1; 
    %     JAC(0,1) = 1;
    %     JAC(1,1) = parameters[1]^2;
    jacobian_sparse = '[]'
    % if max_ord is zero jacobian_sparse will equal '[]'
    % if output_type is anything other than 'cvode_sparse'
    %     jacobian_sparse will equal '[]'
    % 
    % if output_type is equal to 'cvode_sparse'
    %   see also documentation on SUNMatrix_sparse in the CVODES manual
    %   the matrix storage format is compressed sparse row
    %   jacobian_sparse for the example input will be:
    %   'data[0] = parameters[0]^2;
    %    data[1] = 1;
    %    data[2] = 1;
    %    data[3] = parameters[1]^2;
    %    indexvals[0] = 0;
    %    indexvals[1] = 1;
    %    indexvals[2] = 0;
    %    indexvals[3] = 1;
    %    indexptrs[0] = 0;
    %    indexptrs[1] = 2;
    %    indexptrs[2] = 4;'
    
    jacobian_nnz                      = '[]'
    % number of nonzeros of the Jacobian of the right hand side of the system
    % w.r.t. the variables.

    
    jacobian_handle                   = '[]'
    jacobian_params                   = '[]'
    jacobian_params_handle            = '[]'
    d_sensitivity_dt_code     
    d_sensitivity_dt_pars_code
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
    
    analytic_jacobian                 = 'false'
    jacobian_storage                  = 'DENSE'
    jacobian_lower_bandwidth          = 'NOT_APPLICABLE'
    jacobian_upper_bandwidth          = 'NOT_APPLICABLE'
  end

  methods
    % constructor
    %
    %% Inputs:
    % 
    % name: a char array or string containing a filename for the generated
    % system
    % 
    % parameters_str: a char array or string containing the names of the
    % parameters of the system of the ODEs, separated by any number of spaces
    % and/or commas.
    %
    % time: a char array or string containing the name of the variable that
    % represents time
    %
    % max_ord_derivatives: the maximum order of derivatives for which code is
    % generated that evaluates the derivatives
    %
    % equations: a 1d cell array of character vectors containing at least 1
    % differential equation and optionally partial expressions for example:
    % 
    % f = x + y  (optional partial expression)
    % x' = f * y (differential equation)
    % y' = f * x (differential equation)
    %
    % app: an object that implements the method updatStatus(status) or an empty
    % array.
    %
    % output_type: must be one of 'odefile', 'C', or 'cvode'; 
    function s = SystemFileGenerator(...
        name,...
        parameters_str,...
        time,...
        max_ord_derivatives,...
        equations,...
        app,...
        output_type)
            
      % set defining properties
      s.name                = strtrim(char(name));
      s.input_pars_str      = strtrim(char(parameters_str));
      s.equations           = strtrim(cellstr(equations));
      s.max_ord_derivatives = max_ord_derivatives;
      s.time                = strtrim(char(time));
      s.app                 = app;
      s.output_type         = output_type;
      
      % set derived properties
      if isempty(s.input_pars_str)
        SystemFileGenerator.throwException( ...
                'You must specify at least one parameter.');
      end
      
      if isempty(equations)
        SystemFileGenerator.throwException( ...
                'You must specify at least one equation.');
      end
      
      s.input_pars = regexp(s.input_pars_str, '( |,)+','split');
      s.parse_equations();
      
      s.parameter_arguments = s.generate_parameter_arguments();
      s.input_pars_str = strjoin(s.input_pars, ' ');
      s.input_vars_str = strjoin(s.input_vars, ' ');

      
      s.verify_inputs()
      
      syms(s.internal_syms{:})
      for i = 1 : length(s.intermediate_expressions)
        evalc(s.intermediate_expressions{i});
      end
      
      s.rhs_with_substitutions = cell(length(s.input_vars), 1);
      
      for i = 1 : length(s.rhs)
        s.rhs_with_substitutions{i} = char(eval(s.rhs{i}));
      end
        
      s.formatted_rhs = s.format_rhs();
      
      if (s.max_ord_derivatives >= 1)
        s.analytic_jacobian = 'true';
        [jac, storage, uppr_bw, lwr_bw] = s.compute_jacobian('vars');
        s.jacobian                 = jac;
        s.jacobian_storage         = storage;
        s.jacobian_upper_bandwidth = uppr_bw;
        s.jacobian_lower_bandwidth = lwr_bw;
        s.jacobian_params          = s.compute_jacobian('pars');
        % s.compute_sensitivity_right_hand_sides must happen after
        % s.compute_jacobian('vars') and s.compute_jacobian('pars');
        s.compute_sensitivity_right_hand_sides;
        switch s.output_type
          case {'C', 'C_sparse', 'cvode', 'cvode_sparse'}
            s.jacobian_handle       =sprintf('@%s.jacobian_mex'       , s.name);
            s.jacobian_params_handle=sprintf('@%s.jacobian_params_mex', s.name);
          case 'odefile'
            s.jacobian_handle         = '@jacobian';
            s.jacobian_params_handle  = '@jacobian_params';
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
    
    function parse_equations(s)
      left_hand_sides  = cell(length(s.equations), 1);
      right_hand_sides = cell(length(s.equations), 1);


      is_a_DE = false(length(s.equations), 1);
      variables = cell(length(s.equations), 1);
      

      for i = 1 : length(s.input_pars)
        if isempty(regexp(s.input_pars{i},'^[a-zA-Z][a-zA-Z0-9_]*$', 'once'))
          error('Parameter %d ( ''%s'' ) is not a valid Matlab identifier', ...
                  i, s.input_pars{i});
        end
      end

      for i = 1 : length(s.equations)
        parts = strsplit(s.equations{i}, '=');
        if length(parts) ~= 2
          error('SystemFileGenerator:too_much_equal_signs', ...
                  'Equation %d (%s) contains %d equals signs', i, ...
                  s.equations{i}, length(parts) - 1)
        end
        left_hand_sides{i}  = strtrim(parts{1});
        right_hand_sides{i} = strtrim(parts{2});

        invalid_lhs_frmt_strng = [ ...
                  'Equation %d ( %s ) has an invalid left hand side. ' ...
                  'The left hand side must be a valid Matlab identifier, ' ...
                  'optionally followed by an aposthrophe to specify ', ...
                  'a differential equation.' ];

        % if the left hand side ends with a ', then it is a differential
        % equation.
        if endsWith(left_hand_sides{i}, '''')
          is_a_DE(i) = true;
          if isempty(regexp(left_hand_sides{i}, ...
                  '^[a-zA-Z][a-zA-Z0-9_]*''$', 'once'))
            % ^ means start of the string
            % [a-zA-Z] means extactly one letter
            % [a-zA-Z0-9_]* means one or more letters, numbers or underscores
            % '' means one apostrophe
            % $ mean end of the string
            error('SystemFileGenerator:invalid_lhs', ...
                    invalid_lhs_frmt_strng, i, s.equations{i});
          end
          variables{i} = left_hand_sides{i}(1:end-1);  
        else
          if isempty(regexp(left_hand_sides{i}, ...
                  '^[a-zA-Z][a-zA-Z0-9_]*$', 'once'))
            error('SystemFileGenerator:invalid_lhs', ...
                    invalid_lhs_frmt_strng, i, s.equations{i});
          end
          variables{i} = left_hand_sides{i};
        end
      end
      s.input_vars               = variables  (  is_a_DE);
      s.input_vars_str           = strjoin(s.input_vars);
      s.input_intermediates      = variables  (~ is_a_DE);
      s.intermediate_expressions = s.equations(~ is_a_DE);

      s.rhs = right_hand_sides(is_a_DE);
      
      [s.input_syms, s.internal_syms, s.output_syms] = s.generate_sym_lists;
      
      s.internal_vars_str = s.input_to_internal(s.input_vars_str);
      s.internal_vars = strsplit(s.internal_vars_str);

      s.internal_pars_str = s.input_to_internal(s.input_pars_str);
      s.internal_pars = strsplit(s.internal_pars_str);      
      
      s.syms_arg = [s.internal_vars_str ' ' s.internal_pars_str];
      
      for i = 1:length(s.rhs)
        s.rhs{i} = s.input_to_internal(s.rhs{i});
      end
      
      for i = 1:length(s.intermediate_expressions)
        s.intermediate_expressions{i} = ...
                s.input_to_internal(s.intermediate_expressions{i});
      end

      
      eval(['syms ' s.syms_arg]);

      error_id = 'SystemFileGenerator:expression_problem';
      format_str = 'There appears to be a problem with the expression ''%s''.';

      

      for i = 1 : length(s.intermediate_expressions)
        try
          evalc(s.intermediate_expressions{i});
        catch my_exception
          switch my_exception.identifier
            case'MATLAB:UndefinedFunction'
              SystemFileGenerator.throw_undefined_var_error(my_exception)
            otherwise
              error(error_id, format_str, s.intermediate_expressions{i});
          end
        end
      end

      for i = 1 : length(s.rhs)
        try
          evalc(s.rhs{i});
        catch my_exception
          switch my_exception.identifier
            case'MATLAB:UndefinedFunction'
              SystemFileGenerator.throw_undefined_var_error(my_exception)
            otherwise
              expression = [s.input_vars{i} ''' = ' s.rhs{i}];
              error(error_id, format_str, expression);
          end
        end
      end
    end
    
    function out = input_to_internal(s, in)
      out = replace_symbols(in, s.input_syms, s.internal_syms);
    end
    
    function generate_file(s)
      switch s.output_type
        case {'C', 'C_sparse', 'cvode', 'cvode_sparse'}
          s.generate_c_files();
        case 'odefile'
          path = SystemFileGenerator.get_cl_matcontL_path();
          templates_dir = fullfile(path, 'SystemFileGenerator', 'templates');
          filename      = fullfile(path, 'Systems', [s.name,'.m']);
          content       = emat2(fullfile(templates_dir, 'system.m.emat'));
          SystemFileGenerator.print_to_file(filename, content);
      end
    end
    
    function generate_c_files(s)
      mc_path        = get_cl_matcontL_path();
      systems_path   = fullfile(mc_path, 'Systems');
      system_path    = fullfile(systems_path, s.name);
      templates_path = fullfile(mc_path, 'SystemFileGenerator', 'templates');
      
      if ~ exist(system_path, 'dir')
        mkdir(system_path)
      end
      
      filename = fullfile(systems_path, [s.name '_mex.m']);
      content = emat2(fullfile(templates_path, 'system_mex.m.emat'));
      SystemFileGenerator.print_to_file(filename, content);
      
      filename = fullfile(system_path, 'dydt_mex.c');
      content = emat2(fullfile(templates_path, 'system_dydt.c.emat'));
      SystemFileGenerator.print_to_file(filename, content);
      
      if s.max_ord_derivatives >= 1  
        if endsWith(s.output_type, '_sparse')
          content = emat2(fullfile(templates_path, ...
                              'system_jacobian_sparse.c.emat'));
        else
          content = emat2(fullfile(templates_path, 'system_jacobian.c.emat'));
        end

        filename = fullfile(system_path, 'jacobian_mex.c');
        SystemFileGenerator.print_to_file(filename, content);

        filename = fullfile(system_path, 'jacobian_params_mex.c');
        content = emat2(fullfile(templates_path, ...
                          'system_jacobian_params.c.emat'));
        SystemFileGenerator.print_to_file(filename, content);
      end
      
      if any(strcmp({'cvode', 'cvode_sparse'}, s.output_type))

        filename = fullfile(system_path, 'dydt_cvode.c');
        content = emat2(fullfile(templates_path, 'dydt_cvode.c.emat'));
        SystemFileGenerator.print_to_file(filename, content);


        filename = fullfile(system_path , 'user_data.h');
        content = emat2(fullfile(templates_path, 'user_data.h.emat'));
        SystemFileGenerator.print_to_file(filename, content);

        copyfile(fullfile(templates_path, 'set_precision.h'), ...
                 fullfile(system_path,    'set_precision.h'));

        if s.max_ord_derivatives >= 1
          filename = fullfile(system_path, 'jacobian_cvode.c');
          switch s.output_type
            case 'cvode'
              template = 'jacobian_cvode.c.emat';
            case 'cvode_sparse'
              template = 'jacobian_cvode_sparse.c.emat';
          end
          content = emat2(fullfile(templates_path, template));
          SystemFileGenerator.print_to_file(filename, content);

          filename = fullfile(system_path , 'd_sensitivity_dt.c');
          content = emat2(fullfile(templates_path, 'd_sensitivity_dt.c.emat'));
          SystemFileGenerator.print_to_file(filename, content);

          filename = fullfile(system_path , 'd_sensitivity_dt_pars.c');
          content = emat2(fullfile(templates_path, ...
                            'd_sensitivity_dt_pars.c.emat'));
          SystemFileGenerator.print_to_file(filename, content);
        end
      end
      
      filename = fullfile(system_path, 'compile.m');      
      content = emat2(fullfile(templates_path, 'compile.m.emat'));
      SystemFileGenerator.print_to_file(filename, content);
              
      run(filename);
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
      elseif ~ any(strcmp({'odefile', 'C', 'C_sparse', ...
               'cvode', 'cvode_sparse'},  s.output_type))
        s.throwException([s.output_type ' is not a valid output type']);
      elseif startsWith(s.output_type, 'C') && (s.max_ord_derivatives > 1)
        s.throwException( ...
         'C-output of derivatives of order greater than 1 is not implemented.');
      elseif startsWith(s.output_type, 'cvode') && (s.max_ord_derivatives > 1)
        s.throwException( ...
         ['cvode-output of derivatives of order greater than 1 is not ' ...
          'implemented.']);
     
      elseif floor(s.max_ord_derivatives) ~= s.max_ord_derivatives ...
             || s.max_ord_derivatives < 0 || s.max_ord_derivatives > 5
        s.throwException(['max_ord_derivatives must be a integer' ...
                                'which is at least to 0 and at most 5'])
      end
    end

    function str = generate_parameter_arguments(s)
      new_parameters = strcat('par_', s.input_pars);
      str = strjoin(new_parameters, ', ');
    end

    function formatted_rhs = format_rhs(s)
      fprintf('formatting right hand side:\n');
      clear('print_temp', 'time_print_temp');
      switch s.output_type
        case {'C', 'C_sparse', 'cvode', 'cvode_sparse'}
          eval(['syms ' s.syms_arg]);
          rhs_elements = cell(length(s.rhs), 1);
          for i = 1 : length(s.rhs)
            raw_c_code = eval(sprintf('ccode(%s)',s.rhs{i}));
            print_temp('equation %d of %d', i, length(s.rhs));
            tokens    = regexp(raw_c_code, '.*=(.*);', 'tokens');
            c_code    = tokens{1}{1};
            rhs_elements{i} = sprintf('  dydt[%d] = %s;', i-1, c_code);
          end
          formatted_rhs = strjoin(rhs_elements, '\n');
        case {'odefile'}
          formatted_rhs = '';
          if ~ isempty(s.intermediate_expressions)
            formatted_rhs = replace_symbols(...
                    strjoin(s.intermediate_expressions, ';\n'), ...
                    s.input_syms, s.internal_syms);
            formatted_rhs = [formatted_rhs ';\n'];
          end
          formatted_rhs = sprintf([formatted_rhs, ...
                    'dydt =[\n  ', strjoin(s.rhs, ';\n  ') ']']);
      end
      
      print_temp('replacing symbols');
      formatted_rhs = replace_symbols(formatted_rhs, ...
                              s.internal_syms, s.output_syms);
      clear('time_print_temp');
      print_temp('');
    end

    function [in, internal, out] = generate_sym_lists(s)
      symbols_length = length(s.input_vars) + length(s.input_pars) + ...
                       length(s.input_intermediates);
      in       = cell(1, symbols_length);
      internal = cell(1, symbols_length);
      out      = cell(1, symbols_length);
      for i=1:length(s.input_vars)
        in{i}       = s.input_vars{i};
        internal{i} = ['v_', s.input_vars{i}];
        switch s.output_type
          case {'C', 'C_sparse', 'cvode', 'cvode_sparse'}
          out{i}          = sprintf('y[%d]', i-1);
          case 'odefile'
          out{i}          = sprintf('y(%d)', i  );
        end
      end
      vars_len = length(s.input_vars);
      for i = 1:length(s.input_pars)
        in{vars_len + i}       = s.input_pars{i};
        internal{vars_len + i} = ['p_', s.input_pars{i}];
        switch s.output_type
          case {'C', 'C_sparse', 'cvode', 'cvode_sparse'}
          out{vars_len + i}          = sprintf('parameters[%d]', i - 1);
          case 'odefile'
          out{vars_len + i}          = ['par_', s.input_pars{i}]; 
        end
      end
      offset = length(s.input_vars) + length(s.input_pars);
      for i = 1:length(s.input_intermediates)
        in      {offset + i} = s.input_intermediates{i};
        internal{offset + i} = ['i_', s.input_intermediates{i}];
        out     {offset + i} = ['i_', s.input_intermediates{i}];
      end
    end
    
    function [jac, storage, uppr_bw, lwr_bw] = compute_jacobian(s, vars_or_pars)
      storage = 'DENSE';
      uppr_bw = 'NOT_APPLICABLE';
      lwr_bw = 'NOT_APPLICABLE';
      switch vars_or_pars
        case 'vars'; differentiation_vars = s.internal_vars_str;
        case 'pars'; differentiation_vars = s.internal_pars_str;
      end
      my_rhs = strjoin(s.rhs_with_substitutions, ',');
      eval_str = sprintf('jacobian([%s],[%s])', my_rhs, differentiation_vars);
      eval(['syms ' s.syms_arg]);
      jacobian_sym = eval(eval_str);
      fprintf('simplifying Jacobian matrix\n');
      clear('print_temp', 'time_print_temp');
      for i = 1 : numel(jacobian_sym)
        print_temp('element %d of %d ', i, numel(jacobian_sym))
        if jacobian_sym(i) ~= 0
          jac_element_simplified = simplify(jacobian_sym(i), 'All', true, ...
                                            'Seconds', 1);
          jacobian_sym(i) = jac_element_simplified(end);
        end
      end
      clear('time_print_temp');
      print_temp('');
      switch vars_or_pars
        case 'vars'; s.jacobian_vars_sym = jacobian_sym;
        case 'pars'; s.jacobian_pars_sym = jacobian_sym;
      end
      
      switch s.output_type
        
        case {'C', 'C_sparse', 'cvode', 'cvode_sparse'}
        fprintf('formatting jacobian:\n');
        jac_code = SystemFileGenerator.to_c_code_matrix(jacobian_sym);
        
        if strcmp(vars_or_pars, 'vars') && endsWith(s.output_type, '_sparse')
          storage = 'SPARSE';
          s.assemble_sparse_jacobian(jacobian_sym, jac_code)
        end
        if strcmp(vars_or_pars, 'vars') && strcmp(s.output_type, 'cvode')
          sparsity_pattern = cellfun(@(e) ~ isempty(e), jac_code);
          [lwr_bw, uppr_bw] = bandwidth(double(sparsity_pattern));
          if lwr_bw + uppr_bw < length(s.input_vars) / 2
            storage = 'BANDED';
          end
        end
        lines = cell(numel(jacobian_sym),1);
        j = 1;
        for i = 1 : numel(jacobian_sym)
          if ~ isempty(jac_code{i})
            [row, col] = ind2sub(size(jac_code), i);
            lines{j} = sprintf('  JAC(%d,%d) = %s;', row-1, col-1, jac_code{i});
            j = j + 1;
          end
        end
        number_of_nonzeros = j - 1;
        jac = strjoin(lines(1:number_of_nonzeros), '\n');
        
        case 'odefile'
        jac = SystemFileGenerator.symbolic_mat_to_matlab_code(jacobian_sym);
        
      end % of switch s.output_type
      jac = replace_symbols(jac, s.internal_syms, s.output_syms);
    end
    
    function assemble_sparse_jacobian(s, jacobian_sym, jac_code)
      s.jacobian_nnz = nnz(jacobian_sym);
    
      system_size = length(s.input_vars);
      code      = cell (s.jacobian_nnz, 1);
      indexvals = zeros(s.jacobian_nnz, 1);
      indexptrs = zeros(system_size + 1,  1);
      index = 0;
      for j = 1 : system_size
        indexptrs(j) = index;
        for i = 1:system_size
          if jacobian_sym(i,j) ~= 0
            code(index + 1) = jac_code(i,j);
            indexvals(index+1) = i - 1;
            index = index + 1;
          end
        end
      end
      indexptrs(system_size + 1) = index;

      lines = cell(2 * s.jacobian_nnz + system_size + 1,1);
      for i = 1 : s.jacobian_nnz
        lines{i} = sprintf('  data[%d] = %s;', i - 1, code{i});
      end
      for i = 1 : s.jacobian_nnz
        lines{s.jacobian_nnz + i} = sprintf('  indexvals[%d] = %d;', ...
                i - 1, indexvals(i));
      end
      for i = 1 : system_size + 1
        lines{2 * s.jacobian_nnz + i} = sprintf( ...
                '  indexptrs[%d] = %d;', i - 1, indexptrs(i));
      end
      s.jacobian_sparse = replace_symbols(strjoin(lines, '\n'), ...
                                  s.internal_syms, s.output_syms);
    end
    
    function compute_sensitivity_right_hand_sides(s)
      eval(['syms ' s.syms_arg]);
      n_equations = length(s.input_vars);
      
      sensitivity = sym('s',[n_equations, 1]);      
      d_sensitivity_dt = s.jacobian_vars_sym * sensitivity;
      
      ds_dt_p_code = SystemFileGenerator.to_c_code(d_sensitivity_dt);
      lines = cell(n_equations, 1); 
      for i = 1 : n_equations
        lines{i} = sprintf('  ds[%d] = %s;', i-1, ds_dt_p_code{i});
      end

      s.d_sensitivity_dt_code = strjoin(lines, '\n');
      
      % generate lists of char array representations sensitivity variables
      % internal representations: 's1', 's2', ... ( must be valid Matlab names )
      % output representations: 's[0]', 's[1]', ... ( C array locations )
      sensitivity_internal = cell(1, length(s.input_vars));
      sensitivity_output   = cell(1, length(s.input_vars));
      for i = 1 : n_equations
        sensitivity_internal{i} = sprintf('s%d', i);
        sensitivity_output{i}   = sprintf('s[%d]',i-1);
      end
          
      n_params = length(s.input_pars);
      ds_dt_p_code_strings = cell(n_params, 1);
      for i = 1 : n_params
        ds_dt_pars = d_sensitivity_dt + s.jacobian_pars_sym(:,i);
        ds_dt_p_code = SystemFileGenerator.to_c_code(ds_dt_pars);
        lines = cell(n_equations + 2, 1);
        lines{1} = sprintf('  case %d:', i - 1);
        for j = 1 : n_equations
          lines{j+1} = sprintf('  ds[%d] = %s;', j-1, ds_dt_p_code{j});
        end
        lines{end} = '  break;\n';
        ds_dt_p_code_strings{i} = strjoin(lines, '\n');
      end
      
      s.d_sensitivity_dt_pars_code = sprintf([ ...
        '  switch (parameter_index) {\n', ...
        strjoin(ds_dt_p_code_strings, ''), ...
        '  }\n', ...
      ]);

%     [cell_array_1(:); cell_array_2(:)] concatenates cell_array_1 and
%     cell_array_2 into a column cell vector.
      my_internal_symbols = [s.internal_syms(:); sensitivity_internal(:)];
      my_output_symbols   = [s.output_syms(:)  ; sensitivity_output(:)  ];


      s.d_sensitivity_dt_code = replace_symbols(s.d_sensitivity_dt_code, ...
              my_internal_symbols, my_output_symbols);
                    
      s.d_sensitivity_dt_pars_code = replace_symbols( ...
              s.d_sensitivity_dt_pars_code, ...
              my_internal_symbols, my_output_symbols);
    end

    function d = compute_hessians(s, vars_or_pars)
      eval(['syms ' s.syms_arg]);
      switch vars_or_pars
        case 'vars'; diff_vars = s.internal_vars;
        case 'pars'; diff_vars = s.internal_pars;
      end
      len   = length(diff_vars);
      d_sym = sym(zeros(length(s.rhs),len,len));
      for i = 1:length(s.rhs_with_substitutions)
        for j = 1:length(diff_vars)
          for k = 1:j
            evalstr = sprintf('diff(%s, %s, %s)', ...
                    s.rhs_with_substitutions{i}, diff_vars{j}, diff_vars{k});
            pd = eval(evalstr);
            d_sym(i,j,k) = pd;
            d_sym(i,k,j) = pd;
          end
        end
      end
      d = SystemFileGenerator.symbolic_mat_to_matlab_code(d_sym);
      d = replace_symbols(d, s.internal_syms, s.output_syms);
    end
    
    function d = compute_3rd_ord_derivatives(s)
      eval(['syms ' s.syms_arg]);
      variables = s.internal_vars;
      len   = length(variables);
      d_sym = sym(zeros(length(s.rhs),len,len,len));
      for i = 1:length(s.rhs_with_substitutions)
        for j = 1:length(variables)
          for k = 1:j
            for l = 1:k
              evalstr = sprintf('diff(%s, %s, %s, %s)', ...
                  s.rhs_with_substitutions{i}, variables{j}, ...
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
      for i = 1:length(s.rhs_with_substitutions)
        for j = 1:length(variables)
          for k = 1:j
            for l = 1:k
              for m = 1:l
                evalstr = sprintf('diff(%s, %s, %s, %s, %s)', ...
                    s.rhs_with_substitutions{i}, ...
                    variables{j}, variables{k}, variables{l}, ...
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
      d = SystemFileGenerator.symbolic_mat_to_matlab_code(d_sym);
      d = replace_symbols(d, s.internal_syms, s.output_syms);
    end

    function d = compute_5th_ord_derivatives(s)
      eval(['syms ' s.syms_arg]);
      variables = s.internal_vars;
      len = length(variables);
      d_sym = sym(zeros(length(s.rhs),len,len,len,len,len));
      for i = 1:length(s.rhs_with_substitutions)
        for j = 1:length(variables)
          evalstr = sprintf('diff(%s, %s)', s.rhs_with_substitutions{i}, ...
                            variables{j});
          % pd1, pd2, pd3 will be used via the "eval" function
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
      d = SystemFileGenerator.symbolic_mat_to_matlab_code(d_sym);
      d = replace_symbols(d, s.internal_syms, s.output_syms);
    end
    
    function show_status(s,status)
      if ~isempty(s.app)
        s.app.updateStatus(status);
      else
        disp(status)
      end
    end
     
  end
  
  methods (Access = protected)
    
    % Generates a pretty string representation of a SystemFileGenerator object.
    % This string representation is displayed when disp is called with a
    % SystemFileGenerator as argument, or when a SystemFileGenerator object is
    % created in a statement without a semicolon.
    function displayScalarObject(s)
      lines{1} = sprintf('%30s: %s', 'name', s.name);
      lines{2} = sprintf('%30s: %s', 'variables', s.input_vars_str);
      lines{3} = sprintf('%30s: %s', 'parameters', s.input_pars_str);
      lines{4} = sprintf('%30s: %d', 'maximum order of derivatives', ...
                                                         s.max_ord_derivatives);
      lines{5} = sprintf('%30s: %s', 'time variable', s.time);
      lines{6} = sprintf('%30s: %s', 'equations', s.equations{1});
      for i = 2 : length(s.equations)
        lines{6 - 1 + i} = sprintf('%30s  %s', '', s.equations{i});
      end
      lines{end + 1} = sprintf('%30s: %s', 'output type', s.output_type);
      
      truncate = @SystemFileGenerator.truncate;
      if s.max_ord_derivatives >= 1
        lines{end + 1} = sprintf('%30s: %s', 'jacobian', truncate(s.jacobian));
        lines{end + 1} = sprintf('%30s: %s', 'jacobian w.r.t. parameters', ...
                                 truncate(s.jacobian_params));
      end
      if s.max_ord_derivatives >= 2
        lines{end + 1} = sprintf('%30s: %s', 'hessians', truncate(s.hessians));
        lines{end + 1} = sprintf('%30s: %s', 'hessians w.r.t. parameters', ...
                                 truncate(s.hessians_params)); 
      end
      if s.max_ord_derivatives >= 3
        lines{end + 1} = sprintf('%30s: %s', '3rd order derivatives', ...
                                 truncate(s.third_order_derivatives));
      end
      if s.max_ord_derivatives >= 4
        lines{end + 1} = sprintf('%30s: %s', '4th order derivatives', ...
                                 truncate(s.fourth_order_derivatives));
      end
      if s.max_ord_derivatives >= 5
        lines{end + 1} = sprintf('%30s: %s', '5th order derivatives', ...
                                 truncate(s.fifth_order_derivatives));
      end
      out = strjoin(lines, '\n');
      fprintf(out);
      fprintf('\n');
    end
    
    function out = intermediate_to_input(s, in)
      out = replace_symbols(in, s.intermediate_syms, s.input_syms);
    end

  end
    
  methods(Static)
    
    function out = truncate(in)
      truncate_at = 1000;
      if length(in) >= truncate_at
        out = [in(1:truncate_at) ' ... Output was truncated'];
      else
        out = in;
      end
    end
    
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
      exception = MException('SystemFileGenerator:BadInput', msg);
      throw(exception);
    end
    
    function throw_undefined_var_error(UndefinedFunctionError)
      regex = '''\w*''';
      matches = regexp(UndefinedFunctionError.message, regex, 'match');
      variable = matches{1};
      message = sprintf([ ...
                'The system equations contain the variable %s ' ...
                'which you did not specify as a variable in the system ' ...
                'or as a parameter. ' ...
                'Define the time-derivative of %s, declare %s ' ...
                'as a parameter, or create an intermediate expression with ' ...
                'that defines %s and try again.'], ...
                variable, variable, variable, variable);
      undefined_func_error_id = 'SystemFileGenerator:UndefinedFunctionInEval';
      MException(undefined_func_error_id, message).throw;
    end
    
    function s = new(name, par_str, time, max_ord, equations, output_type)
      if nargin == 5
        output_type = 'odefile';
      end
      s = SystemFileGenerator(name, par_str, time, max_ord, ...
              equations, [], output_type);
    end
    
    function print_to_file(filename, content)
      file_id = fopen(filename, 'w');
      fprintf(file_id, '%s', content);
      fclose(file_id);
    end
    
    % input type:  expressions : array of symbolic expressions of any size.
    % output type: symbolic_mat: 1d char array.
    %
    % converts an array of symbolic expressions whose symbols are listed in the
    % instance variable "syms_arg" to a 1D cell array of char arrays containing
    % the C code that evaluates the array of symbolic expressions. Note that the
    % output is 1D even if the input is 2D, 3D or any higher dimension.    
    function c_code = to_c_code(expressions)
      clear('print_temp', 'time_print_temp')
      c_code = cell(numel(expressions),1);
      for i = 1 : numel(expressions)
        if ~ strcmp(char(expressions(i)),'0')
          print_temp('element %d of %d', i, numel(expressions));
          my_c_code = ccode(expressions(i));
          tokens    = regexp(my_c_code, '.*=(.*);', 'tokens');
          c_code{i} = tokens{1}{1};
        end
      end
      clear('time_print_temp')
      print_temp()
    end
        
    % input type:  expressions : array of symbolic expressions of any size.
    % output type: symbolic_mat: 1d char array.
    %
    % converts an array of symbolic expressions whose symbols are listed in the
    % instance variable "syms_arg" to a 1D cell array of char arrays containing
    % the C code that evaluates the array of symbolic expressions. Note that the
    % output is 1D even if the input is 2D, 3D or any higher dimension.    
    function c_code = to_c_code_matrix(expressions)
      clear('print_temp');
      c_code = cell(size(expressions));
      for i = 1 : numel(expressions)
        if ~ strcmp(char(expressions(i)), '0')
          print_temp('element %d of %d', i, numel(expressions));
          my_c_code = ccode(expressions(i));
          tokens    = regexp(my_c_code, '.*=(.*);', 'tokens');
          c_code{i} = tokens{1}{1};
        end
      end
      clear('time_print_temp')
      print_temp();
    end
    
    % input type:  symbolic_mat: array of symbolic expressions of any size.
    % output type: symbolic_mat: 1d char array.
    %
    % converts an array of symbolic expressions whose symbols are listed in
    % the instance variable "syms_arg" to a char array containing
    % the Matlab code that evaluates the array of symbolic expressions.
    %
    % The input "symbolic_mat" should be an array of symbols. "symbolic_mat" can
    % have any dimension.
    function symbolic_mat = symbolic_mat_to_matlab_code(symbolic_mat)
      symbolic_mat = char(matlabFunction(symbolic_mat));
      % The variable "symbolic_mat" now contains a char array with the matlab
      % code that evaluates the input "symbolic_mat". For instance if the input
      % symbolic_mat is [p_a ; p_b], then the variable "symbolic_mat" now
      % contains the string "@(p_a, p_b) [p_a, p_b]".
      
      % We discard everything between the first pair of parentheses by
      % @\(.*?\) and capture the rest by (.*). The question mark makes the
      % expression .* in @\(.*?\) lazy (see mathworks documentation on regexp)
      tokens   = regexp(symbolic_mat, '@\(.*?\)(.*)','tokens');
      symbolic_mat = tokens{1}{1};
    end
      
  end   
  
end
