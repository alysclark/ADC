This folder contains all the .exelem and .exnode files that are outputted by OpenCMISS when running a time dependent problem.
In a time dependent problem, we set the start time of our problem (typically 0) and the end time (i.e. how long we are solving for).
For example, we would set a very long end time in order to reach steady state, or a very small end time in order to see how the concentration gradient changes depending on our initial conditions and boundary conditions.
The size of our time step determines the 'accuracy' of our solution. Typically use a time step that is at least 20x smaller than the end time you are testing.
Any larger and you may notice strange fluctuations when you visualise on cmgui.
The .exnode files contain information on the geometric field (coordinates of each node) and dependent field (concentration at each node). 
The perl file in the ADC->Extra code folder lets you find the average concentration across all the nodes for several time dependent problems, i.e. compare average concentration across mesh at 0.002 seconds and that at 0.004.
This lets you determine whether you have reached/approached steady state.
To run the perl file, just type perl FindAverage.pl
Make sure you have all the files you want to display the average of in the same file as the perl file.
