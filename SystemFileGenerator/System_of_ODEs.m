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
    
    output_type       % string. should equal odefile, C, or cvode
    
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
    sensitivity                       = '[]'
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
        output_type, ...
        varargin)
            
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
        s.analytic_jacobian = 'true';
        [s.jacobian, s.sensitivity,s.jacobian_storage, ...
            s.jacobian_upper_bandwidth, s.jacobian_lower_bandwidth] = ...
                                                     s.compute_jacobian('vars');
        s.jacobian_params = s.compute_jacobian('pars');
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
          path = get_path();

          filename = fullfile(path, '../Systems', [s.name,'.m']);
          file_contents = emat2('system.m.emat');
          fileID = fopen(filename,'w');
          fprintf(fileID,'%s',file_contents);
      end
    end
    
    function generate_c_files(s)
      system = fullfile(get_path(), '../Systems', ['+' s.name]);
      if ~ exist(system, 'dir')
        mkdir(system)
      end
      
      filename = fullfile(system, 'odefile_mex.m');
      contents = emat2('system_mex.m.emat');
      fprintf(fopen(filename,'w'),'%s',contents);
      
      filename = fullfile(system, 'dydt_mex.c');
      contents = emat2('system_dydt.c.emat');
      fprintf(fopen(filename,'w'),'%s',contents);
      
      mex_file = fullfile(system, 'dydt_mex');
      mex(filename,'-output', mex_file);
      
      if (s.max_ord_derivatives < 1)
        return
      end
      
      filename = fullfile(system, 'jacobian_mex.c');
      contents = emat2('system_jacobian.c.emat');
      fprintf(fopen(filename,'w'),'%s',contents);
        
      mex_file = fullfile(system, 'jacobian_mex');
      mex(filename,'-output', mex_file);

      filename = fullfile(system, 'jacobian_params_mex.c');
      contents = emat2('system_jacobian_params.c.emat');
      fprintf(fopen(filename,'w'), '%s', contents);

      mex_file = fullfile(system, 'jacobian_params_mex');
      mex(filename,'-output', mex_file);

    end
    
    function generate_cvode_files(s)
      system    = fullfile(get_path(), '../Systems', ['+' s.name]);
      templates = fullfile(get_path(), 'cvode');
      if ~ exist(system, 'dir')
        mkdir(system)
      end
      
      filename = fullfile(system ,'dydt_cvode.c');
      contents = emat2(fullfile(templates, 'dydt_cvode.c.emat'));
      fprintf(fopen(filename,'w'),'%s',contents);
      
      filename = fullfile(system , 'user_data.h');
      contents = emat2(fullfile(templates, 'user_data.h.emat'));
      fprintf(fopen(filename,'w'),'%s',contents);
      
      if s.max_ord_derivatives >= 1
        filename = fullfile(system ,'jacobian_cvode.c');
        contents = emat2(fullfile(templates, 'jacobian_cvode.c.emat'));
        fprintf(fopen(filename,'w'),'%s',contents);

        filename = fullfile(system , 'd_sensitivity_dt.c');
        contents = emat2(fullfile(templates, 'd_sensitivity_dt.c.emat'));
        fprintf(fopen(filename,'w'),'%s',contents);
      end
      
      
      base_dir = fullfile(get_path(), '..', '..', 'sundails','builddir','src');
      include_dirs = { ...
        fullfile(base_dir, 'cvodes'), ...
        fullfile(base_dir, 'nvector'), ...
        fullfile(base_dir, 'sunlinsol'), ...
        fullfile(base_dir, 'sunnonlinsol'), ...
        fullfile(base_dir, 'sunmatrix') ...
      };
      includes = {};
      for i = 1 : length(include_dirs)
        files = dir(fullfile(include_dirs{i},'**','*.a'));
        files = arrayfun(@(f) fullfile(f.folder,f.name), files, ...
                                                        'UniformOutput', false);
                                                        
        % [cell_array_1(:); cell_array_2(:)] concatenates cell_array_1 and
        % cell_array_2 into a column cell vector.
        includes = [includes(:); files(:)];
      end
      
      if s.max_ord_derivatives >= 1
        mex_arguments = strjoin([
          {fullfile(templates, 'cvode_mex.c')}, ...
          {fullfile(system, 'dydt_cvode.c')}, ...
          {fullfile(system, 'jacobian_cvode.c')}, ...
          {fullfile(system, 'd_sensitivity_dt.c')}, ...
          {['-I' system]}, ...
          {'-g'}, ...  % -g enables debugging symbols
          includes(:)', ...
          {'-output'}, ...
          {fullfile(system, 'cvode')}, ...
        ], ''', ... \n  ''');
      else
        mex_arguments = strjoin([
          {fullfile(templates, 'cvode_mex.c')}, ...
          {fullfile(system, 'dydt_cvode.c')}, ...
          {['-I' system]}, ...
          {'-g'}, ...  % -g enables debugging symbols
          includes(:)', ...
          {'-output'}, ...
          {fullfile(system, 'cvode')}, ...
        ], ''', ... \n  ''');
      end
    
      mex_build = sprintf('mex( ... \n  ''%s'' ... \n)',mex_arguments);
      filename = fullfile(system, 'recompile_cvode_mex.m');
      fprintf(fopen(filename,'w'),'%s', mex_build);
      
      eval(mex_build);
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
      switch s.output_type
        case {'C', 'cvode'}
          eval(['syms ' s.syms_arg]);
          dydt_elements = cell(length(s.rhs),1);
          for i = 1 : length(s.rhs)

            raw_c_code = eval(sprintf('ccode(%s)',s.rhs{i}));
            tokens    = regexp(raw_c_code, '.*=(.*);', 'tokens');
            c_code    = tokens{1}{1};
            dydt_elements{i} = sprintf('  dydt[%d] = %s;', i-1, c_code);
          end
          formatted_rhs = strjoin(dydt_elements, '\n');
        case 'odefile'
          formatted_rhs = ['[' strjoin(s.rhs, '; ') ']'];
      end
      formatted_rhs = replace_symbols(formatted_rhs, ...
                                                s.internal_syms, s.output_syms);
    end

    function [in, internal, out] = generate_sym_lists(s)
      symbols_length = length(s.input_vars) + length(s.input_pars);
      in       = cell(1, symbols_length + 1);
      internal = cell(1, symbols_length + 1);
      out      = cell(1, symbols_length + 1);
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
      % replace 0.0 by 0 in output
      % in      {end} = '0'; is needed for getPropertyGroup method
      in      {end} = '0';
      internal{end} = '0.0';
      out     {end} = '0';
    end
    
    function [jacobian, sensitivity, storage, ...
           upper_bandwidth, lower_bandwidth] = compute_jacobian(s, vars_or_pars)
      storage = 'DENSE';
      upper_bandwidth = 'NOT_APPLICABLE';
      lower_bandwidth = 'NOT_APPLICABLE';
      switch vars_or_pars
        case 'vars'; differentiation_vars = s.internal_vars_str;
        case 'pars'; differentiation_vars = s.internal_pars_str;
      end
      eval_str = sprintf('jacobian([%s],[%s])', strjoin(s.rhs, ','), ...
                                                          differentiation_vars);
      jacobian_sym = s.safe_eval_w_syms(eval_str);
      switch s.output_type
        
        case {'C', 'cvode'}
          
        c_code_of_jac_elements = s.to_c_code_matrix(jacobian_sym);
        
        if strcmp(vars_or_pars, 'vars')
          sparsity_pattern = cellfun(@(e) ~ isempty(e), c_code_of_jac_elements);
          [lower, upper] = bandwidth(double(sparsity_pattern));
          lower_bandwidth = lower;
          upper_bandwidth = upper;
          if lower + upper < length(s.input_vars) / 2
            storage = 'BANDED';
          end
        end
        
        jacobian_lines         = cell(numel(jacobian_sym),1);
        j = 1;
        for i = 1 : numel(jacobian_sym)
          if ~ isempty(c_code_of_jac_elements{i})
            [row, col] = ind2sub(size(c_code_of_jac_elements), i);
            jacobian_lines{j} = sprintf('  JAC(%d,%d) = %s;', row-1, col-1, ...
                                                     c_code_of_jac_elements{i});
	          j = j + 1;
          end
        end
        number_of_nonzeros = j;
        jacobian = strjoin(jacobian_lines(1:number_of_nonzeros-1), '\n');
        
        
        eval(['syms ' s.syms_arg]);
        sens = sym('s',[length(s.input_vars),1]); 
     
        
        if strcmp(vars_or_pars, 'vars')
        
          sensitivity = jacobian_sym * sens;

          c_code_sensitivity = s.to_c_code(sensitivity);
          sensitivity_lines = cell(numel(sensitivity),1); 
          for i = 1 : numel(c_code_sensitivity)
            sensitivity_lines{i} = sprintf('  ds[%d] = %s;', i-1, ...
                                                      c_code_sensitivity{i});
          end

          sensitivity = strjoin(sensitivity_lines, '\n');
          sens_internal = cell(1,length(s.input_vars));
          sens_output   = cell(1,length(s.input_vars));
          for i = 1 : length(s.input_vars)
            sens_internal{i} = sprintf('s%d', i);
            sens_output{i}   = sprintf('s[%d]',i-1);
          end
          
          
          % [cell_array_1(:); cell_array_2(:)] concatenates cell_array_1 and
          % cell_array_2 into a column cell vector.
          sensitivity = replace_symbols(sensitivity, ...
                          [s.internal_syms(:); sens_internal(:)], ...
                          [s.output_syms(:)  ; sens_output(:)  ]);
        else
          sensitivity = '[]';
        end
        
        case 'odefile'
          
        jacobian = s.symbolic_mat_to_matlab_code(jacobian_sym);
        sensitivity = '[]';
       
      end
      jacobian = replace_symbols(jacobian, s.internal_syms, s.output_syms);
     
      
     
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
    function symbolic_mat = symbolic_mat_to_matlab_code(s, symbolic_mat)
      eval(['syms ' s.syms_arg]);
      symbolic_mat = char(matlabFunction(symbolic_mat));
      % The variable "symbolic_mat" now contains a char array with the matlab
      % code that evaluates the input "symbolic_mat". For instance if the input
      % symbolic_mat is [p_a ; p_b], then the variable "symbolic_mat" now
      % contains the string "@(p_a, p_b) [p_a, p_b]".
      
      % We ignore the everything between the first pair of parentheses by
      % @\(.*?\) and capture the rest by (.*). The question mark makes the
      % expression .* in @\(.*?\) lazy (see mathworks documentation on regexp)
      tokens   = regexp(symbolic_mat, '@\(.*?\)(.*)','tokens');
      symbolic_mat = tokens{1}{1};
    end
    
    % input type:  expressions : array of symbolic expressions of any size.
    % output type: symbolic_mat: 1d char array.
    %
    % converts an array of symbolic expressions whose symbols are listed in the
    % instance variable "syms_arg" to a 1D cell array of char arrays containing
    % the C code that evaluates the array of symbolic expressions. Note that the
    % output is 1D even if the input is 2D, 3D or any higher dimension.    
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
    
        % input type:  expressions : array of symbolic expressions of any size.
    % output type: symbolic_mat: 1d char array.
    %
    % converts an array of symbolic expressions whose symbols are listed in the
    % instance variable "syms_arg" to a 1D cell array of char arrays containing
    % the C code that evaluates the array of symbolic expressions. Note that the
    % output is 1D even if the input is 2D, 3D or any higher dimension.    
    function c_code = to_c_code_matrix(s, expressions)
      eval(['syms ' s.syms_arg]);
      c_code = cell(size(expressions));
      for i = 1 : numel(expressions)
        if ~ strcmp(char(expressions(i)), '0')
          my_c_code = ccode(expressions(i));
          tokens    = regexp(my_c_code, '.*=(.*);', 'tokens');
          c_code{i} = tokens{1}{1};
        end
      end
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
        props.jacobian        = ...
                replace_symbols(s.jacobian,        s.output_syms, s.input_syms);
        props.jacobian_params = ...
                replace_symbols(s.jacobian_params, s.output_syms, s.input_syms);
      end
      if s.max_ord_derivatives >= 2
        props.hessians        = ...
                        replace_symbols(s.hessians,s.output_syms, s.input_syms);
        props.hessians_params = ...
                replace_symbols(s.hessians_params, s.output_syms, s.input_syms);
      end
      if s.max_ord_derivatives >= 3
        props.third_order_derivatives = replace_symbols(...
                        s.third_order_derivatives, s.output_syms, s.input_syms);
      end
      if s.max_ord_derivatives >= 4
        props.fourth_order_derivatives = replace_symbols(...
                       s.fourth_order_derivatives, s.output_syms, s.input_syms);
      end
      if s.max_ord_derivatives >= 5
        props.fifth_order_derivatives = replace_symbols(...
                        s.fifth_order_derivatives, s.output_syms, s.input_syms);
      end

      propgrp = matlab.mixin.util.PropertyGroup(props);
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