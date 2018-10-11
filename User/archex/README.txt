arch1     -- Standalone solver (displacement loading)
test_arch -- Matcont driver
archcurve -- Matcont problem description file

Notes:

 - All the mesh information is passed back and forth through a global.  It
   would be better to pass this through the options list, but that would
   require re-working the options to allow cell arrays (rather than just
   matrices) to be passed to run_continuer.

 - It would be nice to plot the displaced structure at various singular 
   points, but I'm not sure how to extract the displacement vector at 
   those points.
 
 - the BP in test_arch_BP0 is sometimes not located. After running 
   test_arch_BP0 a few times, the BP is eventually located. This behavior 
   seems very strange