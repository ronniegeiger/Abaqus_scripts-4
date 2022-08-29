from part import *
from material import *
from section import *
from assembly import *
from step import *
from interaction import *
from load import *
from mesh import *
from optimization import *
from job import *
from sketch import *
from visualization import *
from connectorBehavior import *

import os
import math
import csv
import multiprocessing
import time
import sys

########################################---VARIABLES---#######################################


###################### Job flags
#Create job:
# 0 = no
# 1 = yes
createJob = 1
#Save input file
# 0 = no
# 1 = yes
saveInputFile = 1
#Perform data check on input file
# 0 = no
# 1 = yes
dataCheck = 0
#Submit job: 
# 0 = no
# 1 = yes
submitJob = 0

###################### Boundary conditions
#Temperature difference (C)
tempApplied = 75.0
#Displacement flags
#These flags inidicate if displacement is applied in a direction
# False = no displacement
# True = displacement as indicated in the variable for displacement
dispXflag = True
dispYflag = False
dispZflag= False
dispX = 0.045
dispY = 0.0
dispZ = 0.0

###################### Step flags and variables
#Step name
stpName = 'Step_name'
#Step basics
# Inclusion of nonlinear effects: 
# 0 = OFF
# 1 = ON
nonlinearGeom = 0
# Total time period for the step
tPeriod = 1
#Increment size
# Max. temperature change
deltaIncrement = 1
# Initial increment size
inicialIncrement = 0.05
# Minimum increment size
minIncrement = 1e-15
#Maximum number of increments
maxNumIncrement = 100000
#Amplitude variation for loading magnitudes during the step
# 0 = Step
# 1 = Ramp
loadVariation = 0
#Step parameters
nonlinearGeomBool = [OFF, ON]
variationMagnitudes = [STEP, RAMP]

###################### CAD variables
#Name of the analysis model (unchangeable)
modelName = 'Model-1'

######################################---LOG FUNCTION---##########################################

#This function is used to make message printing easier and simpler
def plog(str):
    
    #Print to Abaqus message window
    print(str)
    
    #Print to terminal
    print >> sys.__stdout__, str
    
    #Print to file
    pfile.write(str)
    pfile.flush()
        
##############################################################################################
##############################################################################################
##############################################################################################
#####################################---MAIN PROGRAM---#######################################


######################################---IMPORTING INP---#######################################

start0 = time.time()

#Set some names
jobName = 'VOL0160_CNT01GNP01'
#Name of the file to save print messages
print_file = jobName + '.txt'
pfile = open(print_file, "a")
plog('############################ START ############################\n')

#Name of the input file
inpName = jobName+'.inp'

start = time.time()
#Import model from inp
mdb.ModelFromInputFile(
    inputFileName=inpName, 
    name=modelName)

end = time.time()
plog('Importing model from inp time: {}\n'.format(end-start))

##############################################################################################
##############################################################################################
##############################################################################################
#####################################---MAIN PROGRAM---#######################################
####################################---NEW CONDITIONS---######################################|

###################################---MODIFY BOUNDARY CONDITIONS---######################################

#Change BC

#Name for the BC
bcDispName = 'Disp-BC-1'
        
#Check if displacement along the x-direction is applied
if dispXflag:
    
    #Se the displacement BC along the x-direction
    mdb.models[modelName].boundaryConditions[bcDispName].setValues(u1=dispX)

#Check if displacement along the y-direction is applied
if dispYflag:
    
    #Se the displacement BC along the y-direction
    mdb.models[modelName].boundaryConditions[bcDispName].setValues(u2=dispY)

#Check if displacement along the z-direction is applied
if dispZflag:
    
    #Se the displacement BC along the z-direction
    mdb.models[modelName].boundaryConditions[bcDispName].setValues(u3=dispZ)

plog('New boundary conditions applied\n')

###################################---JOB SUBMISSION---######################################

if createJob == 1:
    
    #Create and submit job using Abaqus default values
    mdb.Job(atTime=None, contactPrint=OFF, description='', echoPrint=OFF, 
        explicitPrecision=SINGLE, getMemoryFromAnalysis=True, historyPrint=OFF, 
        memory=90, memoryUnits=PERCENTAGE, model=modelName, modelPrint=OFF, 
        multiprocessingMode=DEFAULT, name=jobName, nodalOutputPrecision=SINGLE, 
        numCpus=int(multiprocessing.cpu_count()*0.8), numDomains=int(multiprocessing.cpu_count()*0.8), numGPUs=1, queue=None, resultsFormat=ODB, scratch=
        '', type=ANALYSIS, userSubroutine='', waitHours=0, waitMinutes=0)
    
    if saveInputFile == 1:
        
        #Save input file
        start = time.time()
        mdb.jobs[jobName].writeInput(consistencyChecking=OFF)
        plog("Input file written: {} secs.\n".format(time.time()-start))
    
    if dataCheck == 1:
        
        #Perform data check
        start = time.time()
        mdb.jobs[jobName].submit(consistencyChecking=OFF, datacheckJob=True)
        mdb.jobs[jobName].waitForCompletion()
        plog("Data check on input: {} secs.\n".format(time.time()-start))

    if submitJob == 1:
        
        #Submit job and save time of submission
        start = time.time()
        mdb.jobs[jobName].submit(consistencyChecking=OFF)
        
        #Make the Python script to wait for the Abaqus job to finish
        #In this way the script can measure the execution time of the Abaqus simulation
        mdb.jobs[jobName].waitForCompletion()
        
        plog("Time for Job execution: {}\n".format(time.time()-start))

end = time.time()        
plog("Time for Abaqus model: {}\n".format(end-start0))
pfile.close()

