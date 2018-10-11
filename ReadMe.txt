CL_MATCONTL release 2016p0

Before staring make sure MATLAB's current working directory is set to the 
main directory of CL_MATCONTL.

First, the file 'init.m' needs to be executed to set the necessary paths.
Once init has completed you are ready to begin running CL_MATCONTL. If a 
particular problem requires additional subdirectories to function, they 
must be added to the init file manually. 

In the directory '/Tutorial/' there are two series of script m-files that will 
execute examples for use in CL_MATCONT. Any of the files ending in the
numeral 0 (i.e. 'testbruss_BP0.m' ) can be run without any prerequisite 
files.  Files ending in other numbers generally require the previous named 
files to be run first.  (i.e 'testbruss_BP1.m' cannot be run before 
'testbruss_BP0.m' since it requires data from the earlier file)
Since the data is saved, prerequisite files need only be executed once.

the lines:

opt = contset(opt,'logfile',          1);
opt = contset(opt,'DiagnosticsLevel', 3);

in the test run files can be modified to suit your needs.

'logfile' can be set to true or false.  If there is no logfile, all the
diagnostics will output to the screen.  If set, the log file will be 
created in the logs directory of the problem with the testrun name combined
with a timestamp.
'DiagnosticsLevel' controls the level of output produced to screen or to 
the log file.
NOTE:  Output to screen can slow performance greatly.  Using a logfile for 
higher levels is recomended.  

DiagnosticsLevel can be set to the following values:

-INF:   No Output:  Disables output to screen or log-file.

   0:   Standard Output:   Information about singularities and general 
        information about the curve will be displayed.  Level 0 information
        is always output to the screen when using a log-file.

   1:   Additional Information:   Eigenvalues near zero, the dimension
        of the invariant subspace and unstable subspace (CIS)

   2:   Advanced information:   Convergence information on singularities.  
        Additional error messages concerning convergence of singularity 
        locator functions

   3:   Maximal information:   Info. concerning eigenvalues at singularity 
        points.

   4:   Developer information:  Custom messages for developement of new 
        code.  (Not intended for users)

   5:   Debugging Information:  function markers and other debugging 
        information.  (Not intended for users) 

For additional information refer to the CL_MATCONTL Manual.


