# tgaFoam
The ```tgaFoam``` solver is designed for test the biomass pyrolysis kinetics. In this solver, the ```CRECK``` and the ```RAC``` multicomponent reaction schemes are implemented. The reaction schemes can be chosen in the ```constant/coalCloud1Properites``` file. The reaction schemes can be modified in the ```constant/reactionProperties``` file.
For postProcessing, there are 2 python files for doing that. The ```postProcessSaveCsv.py``` will save figure and result to a csv file. The ```dimNum.py``` will calculate the pyrolysis number. Before calculating the pyrolysis number, remember to remove the file called weighted_k in the ```0``` folder.
