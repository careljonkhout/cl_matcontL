function run_all_demos_fast_to_slow(start_from)

  close all

  clear('directories');

  directories(1).name = 'bruss_1d';
  directories(1).files = { ...
    'testbruss_BP0', ...
    'testbruss_BP1', ...
    'testbruss_HP0', ...
    'testbruss_HP1', ...
    'testbruss_LP0', ...
    'testbruss_LP1', ...
    'testbruss_U0', ...
    'testbruss_U1', ...
    'hopf_for_cycle'};

  directories(2).name = 'fusion';
  directories(2).files = { ...
    'testfusion_BT0', ...
    'testfusion_BT1', ...
    'testfusion_BT2', ...
    'fusion_equilibrium_integration', ...
  };

  directories(3).name = 'bruss_1d';
  directories(3).files = { ...
    'cycle_extend_oc', ...
    'cycle_from_hopf_oc', ...
    'cycle_integration_oc', ...
  };

  directories(4).name = 'fusion';
  directories(4).files = { ...
    'fusion_cycles', ...
    'extend_fusion_cycles'
  };

  directories(5).name = 'bruss_1d';
  directories(5).files = { ...
    'cycle_extend_ss', ...
    'cycle_from_hopf_ss', ...
    'cycle_integration_ss', ...
  };
    
    
  directories(6).name = 'bruss_1d';
  directories(6).files = { ...
    'cycle_extend_ms', ...
    'cycle_from_hopf_ms', ...
    'cycle_integration_ms', ...
  };

  path_to_this_script = get_path();

  my_pause off;
  
  if nargin == 0
    start_from = 1;
  end

  try
    for directory = directories(start_from:end)
      for i = 1:length(directory.files)
        disp(' ')
        disp(repmat('+',2,80))
        disp(' ')
        file_path = fullfile(directory.name, directory.files{i});
        disp(['     Now running:         ' file_path])
        disp(' ')
        disp(repmat('+',2,80))
        disp(' ')
        run(fullfile(path_to_this_script, directory.name, directory.files{i}))
        % to make sure all the parts of a plot of a demo end up in the right
        % plot, we use drawnow
        drawnow
        pause(0.3)
      end
    end
  catch exception
    my_pause on;
    rethrow(exception)
  end
end