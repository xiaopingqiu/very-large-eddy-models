Inconpressible very large eddy simulation models for OpenFOAM derived from [caelus](https://bitbucket.org/appliedccm/caelus-contributors/src/729b39ba56c8ee99d3feb1929ed5b006a315b466/src/libraries/turbulenceModels/incompressible/VLES/?at=master). 

For users of OpenFOAM-2.3.x, please switch to branch "OpenFAOM-2.3.x", and users of OpenFOAM-4.x shoule switch to branch "OpenFOAM-4.x". Other version are not tested, but for version lower than 3.0.x, branch "OpenFOAM-2.3.x" should work after some small modification. And for version 3.0.x, branch "OpenFOAM-4.x" should work. 

Usage: 
1. git clone the appropriate branch 
2. At the directory of the code, run `wmake libso` 
3. Before running you case, add an item at the end of the "contronDict" : libs ("libincompressibleVLESModels.so"). 

Warning: all the models implemented here are not seriously validated, use them at you own risks. If you find some bugs, please leave a message.   

Reference:
1. Kelly, R., Hickman, A.R., Shi, K., Morris, S.C., Jemcov, A., 2016. Very Large Eddy Simulation of a Transonic Axial Compressor Stage. 52nd AIAA/SAE/ASEE Jt. Propuls. Conf. 1–13. doi:10.2514/6.2016-4745. 

2. Stephens, D.W., Sideroff, C., Jemcov, A., 2015. A two equation VLES turbulence model with near-wall delayed behaviour. 7th Asia-Pacific Int. Symp. Aerosp. Technol. 25 – 27 Novemb. 2015, Cairns 25–27. 
