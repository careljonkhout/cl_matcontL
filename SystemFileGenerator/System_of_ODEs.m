classdef System_of_ODEs < handle
  % 
  % see README for instructions on how to use this tool to generate
  % system files for (cl_)matcont/matcontL
  %
  % defining properties
  % These properties are supplied by the user and fully define the system.
  properties
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
    
    output_type       % string. should equal odefile, C, or cvode
    
    app               % optional: GUI class that calls System_of_ODEs
                      % if app is specified then status notifications
                      % will be passed to app by
                      % calling app.updatStatus(status)
                      
  end
  
  % These properties are computed from the defining propeties when an
  % instance of this class is created.
  properties
    input_vars % cell array of char arrays
    input_pars % cell array of char arrays
    
    % symbols (i.e. variables and parameters) suplied by the user 
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
    
    internal_vars_str % char array
    internal_pars_str % char array
    
    % symbols in the form that is written to the system file
    % cell array of char arrays 
    output_syms
    
    jacobian_vars_sym
    jacobian_pars_sym
    
    
    % all fields below contain char arrays
    syms_arg   
    formatted_rhs
    parameter_arguments
    jacobian                          = '[]'
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
    % rhs:   may be supplied as either:
    %        - a n by 1 string array
    %        - a 1 by n string array
    %        - a 1D cell array of 1D char arrays.
    %        - a n by m character array 
    %        where:
    %        - n is the number of equations
    %        - m is the length of the longest equation
    %
    % variables_str: a char array or string containing the names of n
    % state variables of the system of ODEs, separated by any number of spaces,
    % and/or commas.
    %   
    % parameters_str: a char array or string containing the names of the
    % parameters of the system of the ODEs, separated by any number of spaces
    % and/or commas.
    %
    % time: a char array of string containing the name of the variable that
    % represents time
    %
    % max_ord_derivatives: the maximum order of derivatives for which Matlab
    % code is generated that evaluates the derivatives
    %
    % app: an object that implements the method updatStatus(status) or an empty
    % array.
    %
    % output_type: must be one of 'odefile', 'C', or 'cvode'; 
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
      
      % s.verify_inputs()

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
          case {'c', 'cvode'}
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
    
    function out = input_to_internal(s, in)
      out = replace_symbols(in, s.input_syms, s.internal_syms);
    end
    
    function generate_file(s)
      switch s.output_type
        case 'C'
          s.generate_c_files();
        case 'cvode'
          s.generate_c_files();
          s.generate_cvode_files();
        case 'odefile'
          path = System_of_ODEs.get_cl_matcontL_path();
          templates_dir = fullfile(path, 'SystemFileGenerator', 'templates');
          filename      = fullfile(path, 'Systems', [s.name,'.m']);
          content       = emat2(fullfile(templates_dir, 'system.m.emat'));
          System_of_ODEs.print_to_file(filename, content);
      end
    end
    
    function generate_c_files(s)
      path          = System_of_ODEs.get_cl_matcontL_path();
      system_dir    = fullfile(path, 'Systems', ['+' s.name]);
      templates_dir = fullfile(path, 'SystemFileGenerator', 'templates');
      if ~ exist(system_dir, 'dir')
        mkdir(system_dir)
      end
      
      filename = fullfile(system_dir, 'odefile_mex.m');
      content = emat2(fullfile(templates_dir, 'system_mex.m.emat'));
      System_of_ODEs.print_to_file(filename, content);
      
      filename = fullfile(system_dir, 'dydt_mex.c');
      content = emat2(fullfile(templates_dir, 'system_dydt.c.emat'));
      System_of_ODEs.print_to_file(filename, content);
      
      mex_file = fullfile(system_dir, 'dydt_mex');
      mex(filename, ['-I' templates_dir], '-output', mex_file);
      
      if (s.max_ord_derivatives < 1)
        return
      end
      
      filename = fullfile(system_dir, 'jacobian_mex.c');
      content = emat2(fullfile(templates_dir, 'system_jacobian.c.emat'));
      System_of_ODEs.print_to_file(filename, content);
        
      mex_file = fullfile(system_dir, 'jacobian_mex');
      mex(filename, ['-I' templates_dir],'-output', mex_file);
      

      filename = fullfile(system_dir, 'jacobian_params_mex.c');
      content = emat2(fullfile(templates_dir, 'system_jacobian_params.c.emat'));
      System_of_ODEs.print_to_file(filename, content);

      mex_file = fullfile(system_dir, 'jacobian_params_mex');
      mex(filename, ['-I' templates_dir], '-output', mex_file);

    end
    
    function generate_cvode_files(s)
      path          = System_of_ODEs.get_cl_matcontL_path();
      system_path    = fullfile(path, 'Systems', ['+' s.name]);
      templates_path = fullfile(path, 'SystemFileGenerator', 'templates');
      cvodes_path    = fullfile(path, 'SystemFileGenerator', 'cvodes');
      if ~ exist(system_path, 'dir')
        mkdir(system_path)
      end
      
      filename = fullfile(system_path, 'dydt_cvode.c');
      content = emat2(fullfile(templates_path, 'dydt_cvode.c.emat'));
      System_of_ODEs.print_to_file(filename, content);
      
      
      filename = fullfile(system_path , 'user_data.h');
      content = emat2(fullfile(templates_path, 'user_data.h.emat'));
      System_of_ODEs.print_to_file(filename, content);
      
      copyfile(fullfile(templates_path, 'set_precision.h'), ...
               fullfile(system_path,    'set_precision.h'));
      
      if s.max_ord_derivatives >= 1
        filename = fullfile(system_path ,'jacobian_cvode.c');
        content = emat2(fullfile(templates_path, 'jacobian_cvode.c.emat'));
        System_of_ODEs.print_to_file(filename, content);

        filename = fullfile(system_path , 'd_sensitivity_dt.c');
        content = emat2(fullfile(templates_path, 'd_sensitivity_dt.c.emat'));
        System_of_ODEs.print_to_file(filename, content);
        
        filename = fullfile(system_path , 'd_sensitivity_dt_pars.c');
        content = emat2(fullfile(templates_path, 'd_sensitivity_dt_pars.c.emat'));
        System_of_ODEs.print_to_file(filename, content);
      end
    
      cvodes_sources = dir(fullfile(cvodes_path,'**','*.c'));
      cvodes_sources = arrayfun(@(f) fullfile(f.folder,f.name), ...
                            cvodes_sources, 'UniformOutput', false);
                                      
      cvodes_sources = strrep(cvodes_sources, cvodes_path, '');
      
      cvodes_sources = cellfun(@(f) ...
              {sprintf('fullfile(cvodes_path, ''%s'')', f)}, cvodes_sources);

      if s.max_ord_derivatives >= 1
        mex_arguments = strjoin([ ...
          {'fullfile(path, ''SystemFileGenerator'', ''templates'', ''cvode_mex.c'')'}, ...
          {'fullfile(path, ''Systems'', dirname, ''dydt_cvode.c'')'}, ...    
          {'fullfile(path, ''Systems'', dirname, ''jacobian_cvode.c'')'}, ...
          {'fullfile(path, ''Systems'', dirname, ''d_sensitivity_dt.c'')'}, ...
          {'fullfile(path, ''Systems'', dirname, ''d_sensitivity_dt_pars.c'')'}, ...
          {'[''-I'' fullfile(path, ''SystemFileGenerator'', ''cvodes'', ''include'')]'}, ...
          {'[''-I'' fullfile(path, ''Systems'', dirname)]'}, ...
          {' ... ''-g'' % uncomment to enable debugging symbols'''}, ... 
          {' ... not including preprocessor definitions'}, ... 
          {' ... ''CFLAGS=$CFLAGS -g3'' % uncomment to enable debugging symbols '}, ... 
          {' ... including preprocessor definitions (for gcc)' }
          cvodes_sources(:)', ...
          {'''-output'''}, ...
          {'fullfile(path, ''Systems'', dirname, ''cvode'')'}, ...
        ], ', ... \n  ');
      else
        mex_arguments = strjoin([ ...
          {'fullfile(path, ''SystemFileGenerator'', ''templates'', ''cvode_mex.c'')'}, ...
          {'fullfile(path, ''Systems'', dirname, ''dydt_cvode.c'')'}, ...
          {'[''-I'' fullfile(path, ''SystemFileGenerator'', ''cvodes'', ''include'')]'}, ...
          {'[''-I'' fullfile(path, ''Systems'', dirname)]'}, ...
          {' ... ''-g'' % uncomment to enable debugging symbols'''}, ... 
          {' ... not including preprocessor definitions'}, ... 
          {' ... ''CFLAGS=$CFLAGS -g3'' % uncomment to enable debugging symbols '}, ...
          {' ... inlucluding preprocessor definitions (for gcc)' }
          {strjoin(cvodes_sources, ', ... \n  ')}, ...
          {'''-output'''}, ...
          {'fullfile(path, ''Systems'', dirname, ''cvode'')'}, ...
        ], ', ... \n  ');
      end
    
      mex_build = sprintf([
        'path        = System_of_ODEs.get_cl_matcontL_path;\n', ...'
        'cvodes_path = fullfile(path, ''SystemFileGenerator'', ''cvodes'');\n', ...
        'dirname     = ''+%s'';\n', ...      
        'mex( ... \n  %s ... \n)'], ...
        s.name, mex_arguments);
      filename = fullfile(system_path, 'recompile_cvode_mex.m');
      recompile = sprintf( ...
              'function recompile_cvode_mex()\n %s\n end', mex_build);
      System_of_ODEs.print_to_file(filename, recompile);
      eval(mex_build);
    end
      
    function verify_inputs(s)
      all_equations = strjoin(s.rhs);
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
      elseif strcmp(s.output_type, 'cvode') && (s.max_ord_derivatives > 1)
        s.throwException( ...
         ['cvode-output of derivatives of order greater than 1 is not ' ...
          'implemented.']);
      elseif ~ any(strcmp({'odefile', 'C', 'cvode'},s.output_type))
        s.throwException([s.output_type ' is not a valid output type']);
      elseif floor(s.max_ord_derivatives) ~= s.max_ord_derivatives ...
             || s.max_ord_derivatives < 0 || s.max_ord_derivatives > 5
        s.throwException(['max_ord_derivatives must be a integer' ...
                                'which is at least to 0 and at most 5'])
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
      fprintf('formatting right hand side:\n');
      nBytes = fprintf('loading syms');
      switch s.output_type
        case {'C', 'cvode'}
          
          eval(['syms ' s.syms_arg]);
          dydt_elements = cell(length(s.rhs),1);
          for i = 1 : length(s.rhs)
            % print backspaces to overwrite previous console output
            fprintf(repmat('\b', 1, nBytes))
            nBytes = fprintf('equation %d of %d', i, length(s.rhs));

            raw_c_code = eval(sprintf('ccode(%s)',s.rhs{i}));
            tokens    = regexp(raw_c_code, '.*=(.*);', 'tokens');
            c_code    = tokens{1}{1};
            dydt_elements{i} = sprintf('  dydt[%d] = %s;', i-1, c_code);
           
          end
          formatted_rhs = strjoin(dydt_elements, '\n');
        case 'odefile'
          formatted_rhs = ['[' strjoin(s.rhs, '; ') ']'];
      end
      fprintf(repmat('\b', 1, nBytes))
      nBytes = fprintf('replacing symbols');
      formatted_rhs = replace_symbols(formatted_rhs, ...
                                                s.internal_syms, s.output_syms);
      fprintf(repmat('\b', 1, nBytes))
    end

    function [in, internal, out] = generate_sym_lists(s)
      symbols_length = length(s.input_vars) + length(s.input_pars);
      in       = cell(1, symbols_length);
      internal = cell(1, symbols_length);
      out      = cell(1, symbols_length);
      for i=1:length(s.input_vars)
        in{i}           = s.input_vars{i};
        internal{i} = ['v_', s.input_vars{i}];
        switch s.output_type
          case {'C', 'cvode'}
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
          case {'C', 'cvode'}
          out{vars_len + i}          = sprintf('parameters[%d]', i - 1);
          case 'odefile'
          out{vars_len + i}          = ['par_', s.input_pars{i}]; 
        end
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
      eval_str = sprintf('jacobian([%s],[%s])', strjoin(s.rhs, ','), ...
                                                          differentiation_vars);
      jacobian_sym = s.safe_eval_w_syms(eval_str);
      switch vars_or_pars
        case 'vars'; s.jacobian_vars_sym = jacobian_sym;
        case 'pars'; s.jacobian_pars_sym = jacobian_sym;
      end
      
      switch s.output_type
        
        case {'C', 'cvode'}
        fprintf('formatting jacobian:\n');
        jac_code = System_of_ODEs.to_c_code_matrix(jacobian_sym);
        
        if strcmp(vars_or_pars, 'vars')
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
          
        jac = System_of_ODEs.symbolic_mat_to_matlab_code(jacobian_sym);
       
      end
      jac = replace_symbols(jac, s.internal_syms, s.output_syms);
    end
    
    function compute_sensitivity_right_hand_sides(s)
      eval(['syms ' s.syms_arg]);
      n_equations = length(s.input_vars);
      
      sensitivity = sym('s',[n_equations, 1]);      
      d_sensitivity_dt = s.jacobian_vars_sym * sensitivity;
      
      ds_dt_p_code = System_of_ODEs.to_c_code(d_sensitivity_dt);
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
        ds_dt_p_code = System_of_ODEs.to_c_code(ds_dt_pars);
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
      for i = 1:length(s.rhs)
        for j = 1:length(diff_vars)
          for k = 1:j
            evalstr = sprintf('diff(%s, %s, %s)', s.rhs{i}, ...
                                                    diff_vars{j}, diff_vars{k});
            pd = eval(evalstr);
            d_sym(i,j,k) = pd;
            d_sym(i,k,j) = pd;
          end
        end
      end
      d = System_of_ODEs.symbolic_mat_to_matlab_code(d_sym);
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
      d = System_of_ODEs.symbolic_mat_to_matlab_code(d_sym);
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
      d = System_of_ODEs.symbolic_mat_to_matlab_code(d_sym);
      d = replace_symbols(d, s.internal_syms, s.output_syms);
    end

    function result = safe_eval_w_syms(s, evalstr)
      eval(['syms ' s.syms_arg]);
      result = eval(evalstr);             
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
      exception = error(msg);
    end
    
    function s = new(name, var_str, par_str, time, max_ord, rhs, output_type)
      if nargin == 6
        output_type = 'odefile';
      end
      s = System_of_ODEs(name,var_str,par_str,time,max_ord,rhs,[],output_type);
    end
    
    function path = get_cl_matcontL_path()
      stack         = dbstack('-completenames');
      path_elements = strsplit(stack(1).file, filesep);
      path          = strjoin(path_elements(1:end-2), filesep);
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
      nBytes = fprintf('loading syms');
      %eval(['syms ' s.syms_arg]);
      c_code = cell(numel(expressions),1);
      for i = 1 : numel(expressions)
        if ~ strcmp(char(expressions(i)),'0')
          fprintf(repmat('\b', 1, nBytes));
          nBytes = fprintf('element %d of %d', i, numel(expressions));
          my_c_code = ccode(expressions(i));
          tokens    = regexp(my_c_code, '.*=(.*);', 'tokens');
          c_code{i} = tokens{1}{1};
        end
      end
      fprintf(repmat('\b', 1, nBytes));
    end
        
    % input type:  expressions : array of symbolic expressions of any size.
    % output type: symbolic_mat: 1d char array.
    %
    % converts an array of symbolic expressions whose symbols are listed in the
    % instance variable "syms_arg" to a 1D cell array of char arrays containing
    % the C code that evaluates the array of symbolic expressions. Note that the
    % output is 1D even if the input is 2D, 3D or any higher dimension.    
    function c_code = to_c_code_matrix(expressions)
      n_bytes = 0;
      c_code = cell(size(expressions));
      for i = 1 : numel(expressions)
        if ~ strcmp(char(expressions(i)), '0')
          fprintf(repmat('\b', 1, n_bytes));
          n_bytes = fprintf('element %d of %d', i, numel(expressions));
          my_c_code = ccode(expressions(i));
          tokens    = regexp(my_c_code, '.*=(.*);', 'tokens');
          c_code{i} = tokens{1}{1};
        end
      end
      fprintf(repmat('\b', 1, n_bytes));
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