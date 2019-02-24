# Cellulolytic_Biofilms_Simulation_2d
This is the code that I have written for the simulation of cellulolytic biofilms. 
The code simulates the time evolution of these biological structures. This includes growth, 
formation of colonies, consumption of nutrients and cell death. Besides the above phenomena, 
the code also simulates the stochastic attachment of bacterial cells to the underlying domain.
For introducing this mechanism, a combination of stochastic differential equations and an 
impulse function has been used. The time integration is implemented via explicit numerical schemes. 
As a result of using small time steps for stability, the simulation can take a long time to complete.
Some parts of the code have been parallelized to achieve faster runtime speeds. 2D VERSION
