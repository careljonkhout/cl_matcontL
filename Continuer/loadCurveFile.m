function loadCurveFile(curve)
global cds contopts

cds.curve = curve;

curvehandles = feval(cds.curve);

cds.curve_func             = curvehandles{ 1};
cds.curve_defaultprocessor = curvehandles{ 2};
cds.curve_options          = curvehandles{ 3};
cds.curve_jacobian         = curvehandles{ 4};
cds.curve_hessians         = curvehandles{ 5}; %not used
cds.curve_testf            = curvehandles{ 6};
cds.curve_userf            = curvehandles{ 7};
cds.curve_process          = curvehandles{ 8};
cds.curve_singmat          = curvehandles{ 9};
cds.curve_locate           = curvehandles{10};
cds.curve_init             = curvehandles{11}; %not used - Slot to initializer function
cds.curve_done             = curvehandles{12};
cds.curve_adapt            = curvehandles{13};

% if ~contopts.CIS_UsingCIS      % DV 2018
%     cds.curve_CIS_first_point    = [];
%     cds.curve_CIS_step           = [];
% else
    cds.curve_CIS_first_point    = curvehandles{14};
    cds.curve_CIS_step           = curvehandles{15};
% end
