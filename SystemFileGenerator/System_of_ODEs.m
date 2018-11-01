classdef System_of_ODEs
  % 
  % see README for instructions on how to use this tool to generate
  % system files for (cl_)matcont/matcontL
  %
  % defining properties
  % These properties are supplied by the user and fully define the system.
  properties (SetAccess=immutable)
    name              % name of the system, used to create a filename
                      % char array
    variables_str     % list of variables of the system
                      % coded as a char array separated by spaces
    parameters_str    % list of parameters of the system
                      % coded as a char array separated by spaces
    time              % name of the variable that represents time
                      % char array
    max_ord_derivatives
                      % the maximum order of derivatives that is 
                      % to be computed
                      % integer
    rhs               % an n by 1 cell array of character vectors
    
    app               % optional: GUI class that calls System_of_ODEs
                      % if app is given then status notifications
                      % will be passed to app by
                      % calling app.updatStatus(status)
  end
  
  % These properties are computed from the defining propeties when an
  % instance of this class is created.
  properties (SetAccess=immutable)
    input_variables  % cell array of char arrays
    input_parameters % cell array of char arrays
    
    % symbols (i.e. variables and parameters) supllied by the user 
    % cell array of char arrays
    % the order corresponds to intermediate_symbols and output_symbols
    % i.e. input_symbols{i} corresponds to intermediate_symbols{i}
    % and output_symbols{i}
    input_symbols                 
    
    % symbols used by this class when computing derivatives.
    % cell array of char arrays
    intermediate_symbols
    intermediate_variables % cell array of char arrays
    
    intermediate_variables_str  % char array
    intermediate_parameters_str % char array
    
    % symbols in the form that is written to the system file
    % cell array of char arrays 
    output_symbols
    
    % all fields below contain char arrays
    syms_arg   
    formatted_rhs
    parameter_arguments
end
properties
    jacobian                        = '[]'
    jacobian_handle                 = '[]'
    jacobian_for_parameters         = '[]'
    jacobian_params_handle          = '[]'
    hessians                        = '[]'
    hessians_handle                 = '[]'
    hessians_for_parameters         = '[]'
    hessians_params_handle          = '[]'  
    third_ord_derivatives           = '[]'
    third_ord_derivatives_handle    = '[]' 
    fourth_ord_derivatives          = '[]'
    fourth_ord_derivatives_handle   = '[]'
    fifth_ord_derivatives           = '[]'
    fifth_ord_derivatives_handle    = '[]'
    
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
        app)
            
      % set defining properties
      s.name                = strtrim(char(name));
      s.variables_str       = strtrim(char(variables_str));
      s.parameters_str      = strtrim(char(parameters_str));
      s.rhs                 = strtrim(cellstr(rhs));
      s.max_ord_derivatives = max_ord_derivatives;
      s.time                = strtrim(char(time));
      s.app = app;
      
            
      % set derived properties
      s.input_variables  = regexp(s.variables_str ,'( |,)+','split');    
      s.input_parameters = regexp(s.parameters_str,'( |,)+','split');
      if isempty(s.input_parameters{1})
        s.input_parameters = {};
        s.parameters_str = '';
        s.parameter_arguments = '';
      else
        s.parameter_arguments = s.generate_parameter_arguments();
        s.parameters_str = strjoin(s.input_parameters, ' ');
      end
      s.variables_str  = strjoin(s.input_variables, ' ');
      [s.input_symbols, s.intermediate_symbols, s.output_symbols] = ...
        s.generate_symbols_lists;
      
      s.verify_inputs()

      % generate strings for system file
      for i=1:length(s.rhs)
        s.rhs{i} = replace_symbols(s.rhs{i}, ...
          s.input_symbols, s.intermediate_symbols);
      end
      
      s.intermediate_variables_str = replace_symbols( ...
        s.variables_str, s.input_symbols, s.intermediate_symbols);
      s.intermediate_variables = strsplit(s.intermediate_variables_str);
      
      if ~isempty(s.input_parameters{1})
        s.intermediate_parameters_str = replace_symbols( ...
          s.parameters_str, s.input_symbols, s.intermediate_symbols);
      end      
      
      s.syms_arg = [s.intermediate_variables_str ' ' ... 
        s.intermediate_parameters_str];
        
      s.formatted_rhs       = s.format_rhs();
      
      if (s.max_ord_derivatives >= 1)
        s.jacobian                = s.compute_jacobian_for_variables();
        s.jacobian_for_parameters = s.compute_jacobian_for_parameters();
        s.jacobian_handle         = '@jacobian';
        s.jacobian_params_handle  = '@jacobian_params';
      end
      if (s.max_ord_derivatives >= 2)
        s.show_status('Computing second order derivatives...')
        s.hessians                = s.compute_hessians_for_variables();
        s.hessians_for_parameters = s.compute_hessians_for_parameters();
        s.hessians_handle         = '@hessians';
        s.hessians_params_handle  = '@hessians_params';
      end
      if (s.max_ord_derivatives >= 3)
        s.show_status('Computing third order derivatives...')
        s.third_ord_derivatives = s.compute_3rd_ord_derivatives;
        s.third_ord_derivatives_handle = '@third_ord_derivatives';
      end
      if (s.max_ord_derivatives >= 4)
        s.show_status('Computing fourth order derivatives...')
        s.fourth_ord_derivatives = s.compute_4th_ord_derivatives;
        s.fourth_ord_derivatives_handle = '@fourth_ord_derivatives';
      end
      if (s.max_ord_derivatives >= 5)
        s.show_status('Computing fifth order derivatives...')
        s.fifth_ord_derivatives = s.compute_5th_ord_derivatives;
        s.fifth_ord_derivatives_handle = '@fifth_ord_derivatives';
      end
    end

%     function generate_file(s)
%       % Load template
%       emat = EMat('system.m.emat');
%       emat.errchk = false;
%       fullpath = mfilename('fullpath');
%       path_wo_filename = fullpath(1:end-length('\System_of_ODEs'));
%       filename = fullfile(path_wo_filename, ...
%         '..','Systems', [s.name,'.m']);
%        % Render to a file
%       emat.render(filename);
%     end
    
    function generate_file(s)
      file_contents = emat2();
      fullpath = mfilename('fullpath');
      path_wo_filename = fullpath(1:end-length('\System_of_ODEs'));
      filename = fullfile(path_wo_filename, ...
        '..','Systems', [s.name,'.m']);
      fileID = fopen(filename,'w');
      fprintf(fileID,'%s',file_contents);
    end
    
    function verify_inputs(s)
      all_equations = join(s.rhs);
      if isempty(s.name)
        s.throwException('You must provide a system name.');
      elseif isempty(s.variables_str)
        s.throwException('You must provide at least one variable.');
      elseif isempty(s.rhs)
        s.throwException('You must specify a right hand side.');
      elseif length(s.rhs) ~= length(s.input_variables)
        msg = 'The number of variables must equal';
        msg = [msg ' the number of equations.'];
        s.throwException(msg);
      elseif contains(all_equations, ';')
        s.throwException('The equations may not contain a semicolon.');
      elseif contains(all_equations, '=')
        s.throwException('The equations may not contain an equals sign.');
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
      new_parameters = strcat('par_', s.input_parameters);
      str = strjoin(new_parameters, ', ');
    end

    function formatted_rhs = format_rhs(s)
      formatted_rhs = ['[' strjoin(s.rhs, '; ') ']'];
      formatted_rhs = replace_symbols(formatted_rhs, ...
        s.intermediate_symbols, s.output_symbols);
    end

    function [in, intermediate, out] = generate_symbols_lists(s)
      vars_len = length(s.input_variables);
      symbols_length = vars_len + length(s.input_parameters);
      in           = cell(1,symbols_length);
      intermediate = cell(1,symbols_length);
      out          = cell(1,symbols_length);
      for i=1:length(s.input_variables)
        in{i}           = s.input_variables{i};
        intermediate{i} = ['v_', s.input_variables{i}];
        out{i}          = ['xxxxx(', int2str(i), ')'];
      end
      for i = 1:length(s.input_parameters)
        in{vars_len + i}           = s.input_parameters{i};
        intermediate{vars_len + i} = ['p_', s.input_parameters{i}]; 
        out{vars_len + i}          = ['par_', s.input_parameters{i}]; 
      end
    end

    function jacobian = compute_jacobian_safe(s, variables_str)
      eval_str = sprintf('jacobian([%s],[%s])', ...
          strjoin(s.rhs, ','), variables_str);
      jacobian = s.safe_eval_w_syms(eval_str);
      jacobian = char(jacobian);
      jacobian = jacobian(8:end-1);
      jacobian = replace_symbols(jacobian, ...
        s.intermediate_symbols, s.output_symbols);
    end
        
    % returns the jacobian matrix as a string of matlab code
    % with the i-th variable replaced with xxxxx(i)
    % and a parameter a replaced with par_a
    function jacobian = compute_jacobian_for_variables(s)
      jacobian = s.compute_jacobian_safe(s.intermediate_variables_str);
      jacobian = sprintf('reshape(%s,%d,%d)', ...
        jacobian, length(s.input_variables), length(s.input_variables));
    end

    function jacobian = compute_jacobian_for_parameters(s)
      jacobian = s.compute_jacobian_safe(s.intermediate_parameters_str);
      jacobian = sprintf('reshape(%s,%d,%d)', ...
        jacobian, length(s.input_variables), length(s.input_parameters));
    end

    % returns the hessian matrices as a string of matlab code
    % with the i-th variable replaced with xxxxx(i)
    % and a parameter a replaced with par_a
    % hessians(i,j,k) equals 
    % D_{variable(j),variable(k)} s.derivative(i)
    function hessians = compute_hessians(s, var_str, len_vars)
      hessians = '[';
      for i = 1:length(s.rhs)
        hessian = s.compute_hessian(s.rhs{i},var_str);
        hessians = strcat(hessians,hessian,', ');
      end
      hessians = strip(hessians);     % removes trailing space
      hessians = strip(hessians,',');
      hessians = strcat(hessians,']');
      hessians = replace_symbols(hessians, ...
        s.intermediate_symbols, s.output_symbols);
      hessians = sprintf('reshape(%s,%d,%d,%d)', hessians, ...
        length(s.input_variables), len_vars, len_vars);
    end

    function hessians = compute_hessians_for_variables(s)
      hessians = s.compute_hessians(s.intermediate_variables_str, ...
        length(s.input_variables));
    end

    function hessians = compute_hessians_for_parameters(s)
      hessians = s.compute_hessians(s.intermediate_parameters_str, ...
        length(s.input_parameters));
    end

    
    % returns the hessian of func w.r.t. the 
    % variables listed in variables_str as a string of matlab code
    % func must contain only symbols listed in s.symbols_str
    function hessian = compute_hessian(s, func, variables_str)
      eval(['syms ' s.syms_arg]);
      eval_str = sprintf('hessian(%s,[%s])', func, variables_str);
      hessian = eval(eval_str);
      hessian = char(hessian);
      hessian = hessian(8:end-1);
    end
    
    function d = compute_3rd_ord_derivatives(s)
      eval(['syms ' s.syms_arg]);
      variables = s.intermediate_variables;
      len = length(variables);
      d_sym= cell(length(s.rhs),len,len,len);
      for i = 1:length(s.rhs)
        for j = 1:length(variables)
          for k = 1:j
            for l = 1:k
              evalstr = sprintf('diff(%s, %s, %s, %s)', ...
                  s.rhs{i}, variables{j}, ...
                  variables{k}, variables{l});
              pd = char(eval(evalstr));
              ps = perms([j,k,l]);
              for m =1:length(ps)
                d_sym{i,ps(m,1), ps(m,2), ps(m,3)} = pd;
              end
            end
          end
        end
      end
      d='[';
      for i = 1:length(s.rhs)
        d = strcat(d,'[');
        for j = 1:length(variables)
          d = strcat(d,'[');
          for k = 1:length(variables)
            d = strcat(d,'[');
            for l = 1:length(variables)
              d = [d, d_sym{i,j,k,l}, ',']; %#ok<AGROW>
            end
            d = System_of_ODEs.strip_comma_and_add_bracket(d);
          end
          d = System_of_ODEs.strip_comma_and_add_bracket(d);
        end
        d = System_of_ODEs.strip_comma_and_add_bracket(d);
      end
      d = System_of_ODEs.strip_comma_and_add_bracket(d);
      d = replace_symbols(d, s.intermediate_symbols, s.output_symbols);
      d = strip(d,',');
      d = sprintf('reshape(%s,%d,%d,%d,%d)', d, len, len, len, len);
    end

    function d = compute_4th_ord_derivatives(s)
      eval(['syms ' s.syms_arg]);
      variables = s.intermediate_variables;
      len = length(variables);
      d_sym= cell(length(s.rhs),len,len,len);
      for i = 1:length(s.rhs)
        for j = 1:length(variables)
          for k = 1:j
            for l = 1:k
              for m = 1:l
                evalstr = sprintf('diff(%s, %s, %s, %s, %s)', ...
                    s.rhs{i}, variables{j}, variables{k}, variables{l}, ...
                    variables{m});
                pd = char(eval(evalstr));
                p = perms([j, k, l, m]);
                for pi =1:length(p)
                  d_sym{i,p(pi, 1), p(pi, 2), p(pi, 3), p(pi, 4)} = pd;
                end
              end
            end
          end
        end
      end
      d='[';
      for i = 1:length(s.rhs)
        d = strcat(d,'[');
        for j = 1:length(variables)
          d = strcat(d,'[');
          for k = 1:length(variables)
            d = strcat(d,'[');
            for l = 1:length(variables)
              d = strcat(d,'[');
              for m = 1:length(variables)
                d = [d, d_sym{i,j,k,l,m}, ',']; %#ok<AGROW>
              end
              d = System_of_ODEs.strip_comma_and_add_bracket(d);
            end
            d = System_of_ODEs.strip_comma_and_add_bracket(d);
          end
          d = System_of_ODEs.strip_comma_and_add_bracket(d);
        end
        d = System_of_ODEs.strip_comma_and_add_bracket(d);
      end
      d = System_of_ODEs.strip_comma_and_add_bracket(d);
      d = replace_symbols(d, s.intermediate_symbols, s.output_symbols);
      d = strip(d,',');
      d = sprintf('reshape(%s,%d,%d,%d,%d,%d)', ...
        d, len, len, len, len, len);
    end

    function d = compute_5th_ord_derivatives(s)
      eval(['syms ' s.syms_arg]);
      variables = s.intermediate_variables;
      len = length(variables);
      d_sym= cell(length(s.rhs),len,len,len,len,len);
      for i = 1:length(s.rhs)
        for j = 1:length(variables)
          evalstr = sprintf('diff(%s, %s)', s.rhs{i}, variables{j});
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
                    pd = '0';
                  else
                    evalstr = sprintf('diff(pd4, %s)',variables{n});
                    pd = char(eval(evalstr));
                  end
                  p = perms([j, k, l, m, n]);
                  for pi =1:length(p)
                    d_sym{i,p(pi,1), p(pi,2), p(pi,3),...
                      p(pi,4), p(pi,5)} = pd;
                  end
                end
              end
            end
          end
        end
      end
      d='[';
      for i = 1:length(s.rhs)
        d = strcat(d,'[');
        for j = 1:length(variables)
          d = strcat(d,'[');
          for k = 1:length(variables)
            d = strcat(d,'[');
            for l = 1:length(variables)
              d = strcat(d,'[');
              for m = 1:length(variables)
                d = strcat(d,'[');
                for n = 1:length(variables)
                  d = [d, d_sym{i,j,k,l,m,n}, ',']; %#ok<AGROW>
                end
                d = System_of_ODEs.strip_comma_and_add_bracket(d);
              end
              d = System_of_ODEs.strip_comma_and_add_bracket(d);
            end
            d = System_of_ODEs.strip_comma_and_add_bracket(d);
          end
          d = System_of_ODEs.strip_comma_and_add_bracket(d);
        end
        d = System_of_ODEs.strip_comma_and_add_bracket(d);
      end
      d = System_of_ODEs.strip_comma_and_add_bracket(d);
      d = replace_symbols(d, s.intermediate_symbols, s.output_symbols);
      d = strip(d,',');
      d = sprintf('reshape(%s,%d,%d,%d,%d,%d,%d)', ...
        d, len, len, len, len, len, len);
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
    


    function out_str = strip_comma_and_add_bracket(in_str)
      out_str = strip(in_str,',');
      out_str = [out_str '],'];
    end
        
    function throwException(msg)
      exception = MException('System_of_ODEs:BadInput', msg);
      throw(exception);
    end
    
    function s=new(name, var_str, par_str, time, max_ord, rhs)
      try
        s=System_of_ODEs(name,var_str,par_str,time,max_ord,rhs,[]);
      catch ex
        s=[];
        if   strcmp(ex.identifier,'System_of_ODEs:BadInput')...
          || strcmp(ex.identifier,'System_of_ODEs:EvalError')...
          || strcmp(ex.identifier,'System_of_ODEs:UndefinedFunctionInEval')
            disp('There seems to be an error in your input.');
            disp(ex.message);
        else
          ex.rethrow();
        end
      end
    end
    
  end   
  
end