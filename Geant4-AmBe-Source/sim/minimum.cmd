universe		= vanilla
simulation_path       = /hiskp3/loepke/geant4/ExperimentalSetup/
executable		= $(simulation_path)/build/ExperimentalSetup
output			= $(simulation_path)/sim/logfiles/log.out.$(Process)
error			= $(simulation_path)/sim/logfiles/log.err.$(Process)
datapath        = $(simulation_path)/sim/results
log			    = condor_explora.log
getenv			= true
nice_user		= TRUE
request_memory  = 1000 MB

arguments   = $(simulation_path)/build/run.mac
queue 1 
