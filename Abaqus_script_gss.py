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

import numpy as np
import os
import math
import csv
import multiprocessing
import time
import sys

########################################---VARIABLES---#######################################

###################### Files
#Name of the .csv file with the GS parameters (needs to be in the same work directory as Abaqus)
#csv_gnpFile = 'gnp_data_702'
csv_gnpFile = 'gnp_data.csv'
#Name of the .csv file with the RVE geometry parameters (needs to be in the same work directory as Abaqus)
#csv_geomFile = 'sample_geom_5'
csv_geomFile = 'sample_geom.csv'

###################### Mesh flags and variables
#Number of elements per side of the RVE
elementsPerSide = 20
#EmbeddedMesh/HostMesh ratio
meshRatio = 0.5
#Mesh the model:
# 0 = no
# 1 = yes
meshModel = 1
#Re-mesh model
#This flag is used to decide whether files are generated with GNP indices for
#those inside and outside the RVE
# 0 = do not generate files
# 1 = generate files
reMeshModel = 0

#Analsys type for the element
# 1 = coupled temperature-displacement
# 2 = coupled thermal-electrical-structural
# 3 = thermal-electric
# 4 = only displacement
elemTypeCode = 4

#Element library:
# 0 = Standard
# 1 = Explicit (Explicit only works with coupled temperature-displacement elements)
selectedLibrary = 0
#Geometric order:
# 0 = Linear
# 1 = Cuadratic
selectedGeomOrder = 0

#Element types
elemTemp_Disp = [C3D8T, C3D6T, C3D4T]
elemTemp_Elec_Disp = [Q3D8, Q3D6, Q3D4]
elemTemp_Elec = [DC3D8E, DC3D6E, DC3D4E]
elemTemp_Disp_Quadratic = [C3D20RT, UNKNOWN_WEDGE, C3D10MT]
elem_Disp = [C3D8R, C3D6, C3D4]

###################### Job flags
#Create job:
# 0 = no
# 1 = yes
createJob = 1
#Save input file
# 0 = no
# 1 = yes
saveInputFile = 0
#Perform data check on input file
# 0 = no
# 1 = yes
dataCheck = 0
#Submit job: 
# 0 = no
# 1 = yes
submitJob = 1

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
dispX = 0.5
dispY = 0.0
dispZ = 0.0

###################### Matrix properties
#Define the name of the matrix
matrixMaterial = 'Polypropylene'
#Matrix name
matrixName = 'MATRIX'
#Define the mass density (kg/m3)
matrixDensity = 900
#Define the elastic modulus (GPa)
matrixModulus = 0.59
#Define the Poisson ratio
matrixPoissonR = 0.42
#Define de coefficient of thermal expansion (e-5 1/C)
matrixExpCoeff = 18
#Define the electrical conductivity (S/m)
matrixElecConductivity = 320e-15
#Define the thermal conductivity (W/m*K) - Because of Abaqus needs, but doesn't affect the thermo-mechanical simulation (could be = 1)
matrixThermConductivity = 0.19
#Define the specific heat (J/mol*K) - Because of Abaqus needs, but doesn't affect the thermo-mechanical simulation (could be = 1)
matrixSpecHeat = 75

###################### Nanofiller properties
#Define the name of the filler
fillerMaterial = 'GRAPHENE'
#Define the mass density (kg/m3)
fillerDensity = 2200
#Define the elastic modulus (GPa) - 1000
fillerModulus = 1029
#Define the Poisson ratio - 0.165
fillerPoissonR = 0.149
#Define de coefficient of thermal expansion (e-5 1/C)
fillerExpCoeff = 0.5
#Define the electrical conductivity (S/m)
fillerElecConductivity = 10e7
#Define the thermal conductivity (W/m*K) - Because of Abaqus needs, but doesn't affect the thermo-mechanical simulation (could be = 1)
fillerThermConductivity = 3000
#Define the specific heat (J/mol*K) - Because of Abaqus needs, but doesn't affect the thermo-mechanical simulation (could be = 1)
fillerSpecHeat = 7
#Maximum dimensions
lxMax = 1.0
lyMax = 1.0
tMax = 0.03
margin = 1.01*sqrt(lxMax*lxMax + lyMax*lyMax + tMax*tMax)

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
inicialIncrement = 0.1
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

#Selecting the desired element library
elementLibrary = [STANDARD, EXPLICIT]

###################### CAD variables
#Bigger box name
biggerBoxName = 'BIGGER_BOX'
#Name of the analysis model (unchangeable)
modelName = 'Model-1'
#Sheetsize for the drawing space
sheetSize = 1.0

###################### Constants
#Conversion factor feom radians to degrees
rad2deg = 180/math.pi

#"Zero" for comparing floating point numbers
Zero = 1e-7

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
    
######################################---STRING FUNCTIONS---########################################

#This function generates a string for the filler part
def string_part(filler, i):
    #Return the string in the format FILLER-i
    return '%s-%d' %(filler, i)

#This function generates a string for the filler instance
def string_instance(filler, i):
    #Return the string in the format FILLER-i-1
    return '%s-%d-1' %(filler, i)

#This function generates a string for the node set of GS i
def gs_string_node_set(gs_i):
    #Return the string in the format GS-nodes-gs_i
    return 'GS-%d-NODES' %(gs_i)
    
def ee_string(filler, i):
    return 'EE-%s-%d'%(filler, i)

###################################---ABAQUS FUNCTIONS---#####################################

#Select the element type based on the choice of geometric order
def Select_Elemet_Type(selectedGeomOrder):
    
    #Identify the geometric order
    if selectedGeomOrder == 0:
        if elemTypeCode == 1:
            return elemTemp_Disp
        elif elemTypeCode == 2:
            return elemTemp_Elec_Disp
        elif elemTypeCode == 3:
            return elemTemp_Elec
        elif elemTypeCode == 4:
            return elem_Disp
        else:
            print('WARNING. Invalid element type. Elment type code was set to displacement elements (elemTemp_Disp)')
            return elem_Disp
    elif selectedGeomOrder == 1:
        if elemTypeCode == 1:
            return elemTemp_Disp_Quadratic
    else:
        print('WARNING. Invalid geometric order. Geometric order was set to linear. Elment type code was set to displacement elements (elemTemp_Disp)')
        return elem_Disp

#Create matrix
def Create_Matrix(model, sheetsz, part, P0, Lxyz):

    #Create the matrix geometry
    mdb.models[model].ConstrainedSketch(name='__profile__', sheetSize=sheetsz)
    #Define the two corners of a rectangle on the xy plane
    mdb.models[model].sketches['__profile__'].rectangle(
        point1=(P0[0], P0[1]), point2=(P0[0]+Lxyz[0], P0[1]+Lxyz[1]))
    #Name the part
    mdb.models[model].Part(dimensionality=THREE_D, name=part, type=DEFORMABLE_BODY)
    
    #Use the length along the z-direction as the extrusion depth
    mdb.models[model].parts[part].BaseSolidExtrude(
        depth=Lxyz[2], sketch=mdb.models[model].sketches['__profile__'])
    #Delete the sketch
    del mdb.models[model].sketches['__profile__']

#Create a box bigger than the RVE, which is used to trim the GNPs
def Create_BiggerBox(model, sheetsz, part, P0, Lxyz, margin):
    
    #The length of the bigger box along each direction is the same as the RVE plus the margin
    lengthBiggerBox_x = Lxyz[0]+2.0*margin
    lengthBiggerBox_y = Lxyz[1]+2.0*margin
    lengthBiggerBox_z = Lxyz[2]+2.0*margin

    #Create the geometry for the bigger box
    #Draw the face on the xy plane using two points of a rectangle
    mdb.models[model].ConstrainedSketch(name='__profile__', sheetSize=sheetsz)
    #Note that coordinates of the min point will be (P0[0]-margin, P0[1]-margin, 0)
    mdb.models[model].sketches['__profile__'].rectangle(
        point1=(P0[0]-margin, P0[1]-margin), 
        point2=(P0[0]-margin+lengthBiggerBox_x, P0[1]-margin+lengthBiggerBox_y))
    mdb.models[model].Part(dimensionality=THREE_D, name=part, type=DEFORMABLE_BODY)

    mdb.models[model].parts[part].BaseSolidExtrude(
        depth=lengthBiggerBox_z, sketch=mdb.models[model].sketches['__profile__'])
    del mdb.models[model].sketches['__profile__']

def Create_Matrix_Instance(modelName, matrixName, P0):

    #Instance name of matrix
    str_mat = matrixName + '-1'

    #Create an instance of the matrix
    mdb.models[modelName].rootAssembly.Instance(dependent=ON, name=str_mat, part=mdb.models[modelName].parts[matrixName])

    #Check if the instance needs to bre translated along the z-axis
    if abs(P0[2]) > Zero:

        #z-coordinate of "min point" of RVE is non-zero
        #Current coordinates of min point are (P0[0], P0[1], 0)
        #Move that corner to (P0[0], P0[1], P0[2])
        #i.e: endpoint - starting point = (P0[0], P0[1], P0[2]) - (P0[0], P0[1], 0)
         mdb.models[modelName].rootAssembly.translate(
             instanceList=(str_mat, ),
             vector=(0.0, 0.0, P0[2]))

#Create the hollow box that is used to cut all GS that are partially outside the RVE
def Create_CuttingBox(model, partMatrix, partBox, P0, margin):
    
    #Create an instance of the bigger box
    mdb.models[model].rootAssembly.Instance(dependent=OFF, name=partBox + '-1', 
        part=mdb.models[model].parts[partBox])

    #Current coordinates of min point for partBox are (P0[0]-margin, P0[1]-margin, 0)
    #Move that corner to (P0[0]-margin, P0[1]-margin, P0[2]-margin)
    #i.e: endpoint - starting point = 
    # (P0[0]-margin, P0[1]-margin, 0) - (P0[0]-margin, P0[1]-margin, P0[2]-margin)
    mdb.models[model].rootAssembly.translate(instanceList=(partBox + '-1', ), 
        vector=(0.0, 0.0, P0[2]-margin))
    
    #Create the hollow box by cutting the matrix off of the bigger box
    mdb.models[model].rootAssembly.InstanceFromBooleanCut(
        cuttingInstances=(
            mdb.models[model].rootAssembly.instances[partMatrix + '-1'], ),
        instanceToBeCut=mdb.models[model].rootAssembly.instances[partBox + '-1'],
        name='CUTTER',
        originalInstances=SUPPRESS)

    mdb.models[model].rootAssembly.features[partMatrix + '-1'].resume()

#Create all parts and instances for GSs
def Create_All_GSs(modelName, fillerMaterial, sheetSize, n_gs, P0, corner):
    
    #Variables to check progress
    check_step = 0.1
    frac_thres = 0.1
    
    #Get the starting time
    start_gss = time.time()

    #Iterate over all GSs in the RVE
    for numPart in range(n_gs):
        
        #Generate the part and instance name of the current GNP
        gs_part_str = string_part('GS', numPart)
        gs_inst_str = string_instance('GS', numPart)
        
        #Create a GS part
        Create_GS(modelName, sheetSize, numPart, gs_part_str)
        
        #Generate GNP instance
        mdb.models[modelName].rootAssembly.Instance(
            dependent=ON, name=gs_inst_str,
            part=mdb.models[modelName].parts[gs_part_str])
        
        #Create a set for each vertex of each GS
        GS_Instances_NodeSet(modelName, gs_part_str, gs_inst_str)
        
        #Move every GS to its right place in the RVE and then cut it (if its partially outside)
        Translate_Rotate_and_Cut_GS(modelName, numPart, P0, corner, gs_part_str, gs_inst_str)
        
        #Assign material to GS
        Assign_Section(modelName, fillerMaterial, gs_part_str)
        
        #Calculate the fraction of generated GS parts
        frac = float(numPart+1)/float(n_gs )
        
        #Check if more than the fraction threshold has been generated
        if frac >= frac_thres:
            
            #Update the threshold
            while frac >= frac_thres:
                
                #Increase the threshold by a step
                frac_thres += check_step
                
            #Send message with percentage generated
            plog('   GS parts and instances generated: {} % ({} secs)\n'.format( int( (frac_thres - check_step)*100 ), time.time()-start_gss))
    
#Create GS
def Create_GS(model, sheetsz, gnp_i, gs_part_str):

    #Get the dimensions needed to draw the GS
    half_x = data_gnp[gnp_i][0]/2.0
    half_y = data_gnp[gnp_i][1]/2.0
    thickness = data_gnp[gnp_i][2]
    
    #Set up the sketch
    mdb.models[model].ConstrainedSketch(name='__profile__', sheetSize=sheetsz)
    #Draw the base of the paralellepiped
    mdb.models[model].sketches['__profile__'].rectangle(
        point1=(-half_x, -half_y),
        point2=(half_x, half_y))
    #Name the part
    mdb.models[model].Part(dimensionality=THREE_D, name=gs_part_str, type=DEFORMABLE_BODY)
    #Extrude the base
    mdb.models[model].parts[gs_part_str].BaseSolidExtrude(
    depth=thickness, sketch=mdb.models[model].sketches['__profile__'])
    #Delete sketch
    del mdb.models[model].sketches['__profile__']

#Create material and assign the properties
def Create_Material(model, materialName, massDensity, elasModulus, poissonRatio, expanCoefficient, thermConductivity, specHeat, elecConductivity):
    
    #Assign the material properties required by Abaqus
    mdb.models[model].Material(name=materialName)
    mdb.models[model].materials[materialName].Density(table=((massDensity, ), ))
    mdb.models[model].materials[materialName].Elastic(table=((elasModulus, poissonRatio), ))
    mdb.models[model].materials[materialName].Expansion(table=((expanCoefficient, ), ))
    mdb.models[model].materials[materialName].Conductivity(table=((thermConductivity, ), ))
    mdb.models[model].materials[materialName].SpecificHeat(table=((specHeat, ), ))
    mdb.models[model].materials[materialName].ElectricalConductivity(table=((elecConductivity, ), ))

#Create material section
def Create_Section(model, materialName):
    mdb.models[model].HomogeneousSolidSection(
        material=materialName, 
        name=materialName, 
        thickness=None)

#Assign section to part
def Assign_Section(model, materialName, part_str):

    mdb.models[model].parts[part_str].SectionAssignment(
        offset=0.0, 
        offsetField='', 
        offsetType=MIDDLE_SURFACE, 
        region=Region(cells=mdb.models[model].parts[part_str].cells),
        sectionName=materialName, 
        thicknessAssignment=FROM_SECTION) 

#Translate and rotate the GS
def Translate_Rotate_and_Cut_GS(model, i, P0, corner, gs_part_str, gs_inst_str):
    
    checkOutside = False

    GS_assembly = mdb.models[model].rootAssembly

    #Get all vertices from a graphene sheet
    GS_vertex = GS_assembly.instances[gs_inst_str].vertices

    #thickness = data_gnp[i][2]
    #y_angle = data_gnp[i][3]*(180/math.pi)
    #z_angle = data_gnp[i][4]*(180/math.pi)
    #final_point = data_gnp[i][5:]

    #Translate GS centroid (x0,y0,z0) to the origin
    #I.e: endpoint - starting point = (0,0,0) - (x0,y0,z0)
    GS_assembly.translate(
        instanceList=(gs_inst_str, ),
        vector=(0.0, 0.0, -data_gnp[i][2]/2.0))

    #Rotation in y-axis
    GS_assembly.rotate(
        angle = data_gnp[i][3]*rad2deg,
        axisDirection=(0.0, 1.0, 0.0), 
        axisPoint=(0.0, 0.0, 0.0), 
        instanceList=(gs_inst_str, ))

    #Rotation in z-axis
    GS_assembly.rotate(
        angle = data_gnp[i][4]*rad2deg,
        axisDirection=(0.0, 0.0, 1.0), 
        axisPoint=(0.0, 0.0, 0.0), 
        instanceList=(gs_inst_str, ))

    #Translate GS centroid (currently in the origin) to its new position (x1,y1,z1)
    #I.e: endpoint - starting point = (x1,y1,z1) - (0,0,0)
    GS_assembly.translate(
        instanceList=(gs_inst_str, ),
        vector=data_gnp[i][5:])

    #Cut the GS:
    
    #Iterate over the GNP vertices
    for numVertex in range(8):
    
        #Get coordinates of a current vertex from GNP
        #I THINK THIS CAN BE SIMPLIFIED
        GS_vertex_coord = GS_assembly.getCoordinates(GS_vertex[numVertex])
        
        #Check if the current vertex is outside the RVE
        checkOutside = Is_Point_Outside_RVE(P0, corner, GS_vertex_coord)
        if checkOutside:
            
            #Since the current vertex is outside the RVE, trim it
            GS_assembly.InstanceFromBooleanCut(
                cuttingInstances=(GS_assembly.instances['CUTTER-1'], ), 
                instanceToBeCut=GS_assembly.instances[gs_inst_str], name=gs_part_str,
                originalInstances=SUPPRESS)
            GS_assembly.features['CUTTER-1'].resume()
            
            #Add current GNP number to the list of GNPs that are partially outside the RVE
            indexOutside.append(i)

            #Termineate the for-loop as it is not needed to check if more
            #vertices are outside the RVE
            break

    if not checkOutside:
        indexInside.append(i)    

    #This seems to be needed for the last GNP
    if i == N_GSs - 1:
        GS_assembly.features['CUTTER-1'].suppress()
 
#This function returns true if a point is outside the RVE
def Is_Point_Outside_RVE(P0, corner, point):
    
    #Check if the point coordinates are inside the RVE defined by P0 and corner
    if point[0] < P0[0] or point[0] > corner[0] or point[1] < P0[1] or point[1] > corner[1] or point[2] < P0[2] or point[2] > corner[2]:
        
        #At least one of the conditions was met, so the point is outside the RVE
        return True
    
    #No condition was met, so the point is inside the RVE
    return False
 
#Temperature boundary conditions
def Add_Boundary_Conditions(model, temp):
    
    #Set name of the amplitude
    ampName = 'Amp-1'
    
    #Create a simple amplitude
    mdb.models[model].TabularAmplitude(
        data=((0.0, 0.0), (1.0, 1.0)), 
        name=ampName, 
        smooth=SOLVER_DEFAULT, 
        timeSpan=STEP)

    #Add the type of step depending on the element type, which also determines the type of study
    if elemTypeCode == 1:
        
        #Create themomechanical step
        mdb.models[model].CoupledTempDisplacementStep(
            timePeriod=tPeriod,
            nlgeom=nonlinearGeomBool[nonlinearGeom],
            deltmx=deltaIncrement,
            initialInc=inicialIncrement,
            name=stpName,
            previous='Initial',
            maxNumInc=maxNumIncrement,
            minInc=minIncrement,
            amplitude = variationMagnitudes[loadVariation])
            
    elif elemTypeCode == 2:
        
        #Create Thermoelectromechanical step
        mdb.models[model].CoupledThermalElectricalStructuralStep(
            timePeriod=tPeriod,
            nlgeom=nonlinearGeomBool[nonlinearGeom],
            deltmx=deltaIncrement,
            initialInc=inicialIncrement,
            name=stpName,
            previous='Initial',
            maxNumInc=maxNumIncrement,
            minInc=minIncrement,
            amplitude = variationMagnitudes[loadVariation])
            
    elif elemTypeCode == 3:
        
        #Create Thermoelectric step
        mdb.models[model].CoupledThermalElectricStep(
            timePeriod=tPeriod,
            nlgeom=nonlinearGeomBool[nonlinearGeom],
            deltmx=deltaIncrement,
            initialInc=inicialIncrement,
            name=stpName,
            previous='Initial',
            maxNumInc=maxNumIncrement,
            minInc=minIncrement,
            amplitude = variationMagnitudes[loadVariation])
        
    elif elemTypeCode == 4:
        
        #Name for the BC
        bcDispName = 'BC-Disp'
        
        #Create mechanical step
        mdb.models[model].StaticStep(
            initialInc=inicialIncrement, 
            name=stpName, 
            noStop=OFF, 
            previous='Initial', 
            timeIncrementationMethod=FIXED)
        
        #Generate an array for the reference points to which PBCs will be applied
        pointsPBCs = []
        
        #Append indicated by the flags
        #Check if displacement along the x-direction is applied
        if dispXflag:
            
            #Append reference point to apply displacement along the x-direction
            pointsPBCs.append(mdb.models[model].rootAssembly.sets['RP_X'].referencePoints[0])
        
        #Check if displacement along the y-direction is applied
        if dispYflag:
            
            #Append reference point to apply displacement along the y-direction
            pointsPBCs.append(mdb.models[model].rootAssembly.sets['RP_Y'].referencePoints[0])
        
        #Check if displacement along the z-direction is applied
        if dispZflag:
            
            #Append reference point to apply displacement along the z-direction
            pointsPBCs.append(mdb.models[model].rootAssembly.sets['RP_Z'].referencePoints[0])
         
        #Apply the displacement BC to 
        mdb.models[model].DisplacementBC(
            amplitude=ampName, 
            createStepName=stpName, 
            distributionType=UNIFORM, 
            fieldName='', 
            fixed=OFF, 
            localCsys=None, 
            name=bcDispName, 
            region=Region(referencePoints=pointsPBCs))
        
        #Set the displacement indicated by the flags
        #Check if displacement along the x-direction is applied
        if dispXflag:
            
            #Se the displacement BC along the x-direction
            mdb.models[model].boundaryConditions[bcDispName].setValues(u1=dispX)
        
        #Check if displacement along the y-direction is applied
        if dispYflag:
            
            #Se the displacement BC along the y-direction
            mdb.models[model].boundaryConditions[bcDispName].setValues(u2=dispY)
        
        #Check if displacement along the z-direction is applied
        if dispZflag:
            
            #Se the displacement BC along the z-direction
            mdb.models[model].boundaryConditions[bcDispName].setValues(u3=dispZ)
    
    #Add Temperature boundary condition if temperature is added
    #This happens when the element code is not 4
    if elemTypeCode != 4:
        mdb.models[model].TemperatureBC(
            createStepName='Initial', 
            distributionType=UNIFORM, 
            fieldName='', 
            magnitude=0.0, 
            name='TEMP',
            region=mdb.models[model].rootAssembly.sets['External_Faces'])
            
        mdb.models[model].boundaryConditions['TEMP'].setValuesInStep(
            magnitude=temp, stepName=stpName)
        
        mdb.models[model].boundaryConditions['TEMP'].setValues(amplitude=ampName)

#Embedded elements
def Create_Embedded_Elements(model, matrixName, gs_i, gs_inst_str):
    
    #Add embedded region to GS i
    mdb.models[model].EmbeddedRegion(
        absoluteTolerance=0.0, 
        embeddedRegion=Region(cells=mdb.models[model].rootAssembly.instances[gs_inst_str].cells),
        fractionalTolerance=0.05, 
        hostRegion=Region(cells=mdb.models[model].rootAssembly.instances[matrixName + '-1'].cells),
        name='Element-%d' %(gs_i),
        toleranceMethod=BOTH, 
        weightFactorTolerance=1e-06)

#Embedded elements of the cutted GS
def Create_Embedded_Elements_CuttedGS(model, matrixName, i):
    
    #Add embedded region to GS i, which lies partially outside the RVE
    mdb.models[model].EmbeddedRegion(
        absoluteTolerance=0.0, 
        embeddedRegion=Region(
            cells=mdb.models[model].rootAssembly.instances['GS-%d-2' %(i)].cells.getSequenceFromMask(mask=('[#1 ]', ), )), 
        fractionalTolerance=0.05, 
        hostRegion=Region(
            cells=mdb.models[model].rootAssembly.instances[matrixName + '-1'].cells.getSequenceFromMask(mask=('[#1 ]', ), )), 
        name='Element-%d' %(i), 
        toleranceMethod=BOTH, 
        weightFactorTolerance=1e-06)

#This function generates all meshes
def Generate_Meshes(modelName, matrixName, selectedElementCode, eeMeshSize, numGSs, indexInside, indexOutside):
    
    #Generate the meshes for the GSs
    
    #Index to iterate over the GSs that are partially outside the RVE
    idx = 0
    
    #Number of GSs partially outside the RVE
    gsOut = len(indexOutside)
    #print('gsOut=%d'%(gsOut))
    
    #Iterate over all GSs
    for gs_i in range(numGSs):

        #Booloean to determine if a GS is partially outside the RVE
        partOut = False
        
        #Check if current GS is inside the RVE or not
        #print('idx=%d'%(idx))
        if gsOut != 0 and gs_i == indexOutside[idx]:
            
            #Set to true the flag for GS partially outside the RVE
            partOut = True
            
            #Increase the index and limit it by the total number of GSs
            #partially outside the RVE
            idx = (idx + 1) % gsOut
            #print('    idx=%d'%(idx))
        
        #Get the part name
        gs_part_str = string_part('GS', gs_i)
        
        #Get the instance name depending on where is the GS
        if partOut:
            
            #Get the instance for a GS that is trimmed
            gs_inst_str = gs_part_str + '-2'
            
        else:
            #GS is inside the RVE
            #Get the instance name using the default function
            gs_inst_str = string_instance('GS', gs_i)
            
        #Create embedded elments constraint
        Create_Embedded_Elements(modelName, matrixName, gs_i, gs_inst_str)
        
        #Mesh each graphene sheet if the flag for meshing is set to 1
        #(meshing: GS mesh element size < Matrix mesh element size)
        if meshModel == 1:
            Mesh_GS(modelName, selectedElementCode, eeMeshSize, gs_part_str)
            
            #Check if:
            #the GS is partially outside the RVE
            #AND
            #(if so) the GS was actually meshed, which is done by checking the number of meshed regions
            if partOut and mdb.models[modelName].parts[gs_part_str].getMeshStats((mdb.models[modelName].parts[gs_part_str].cells,)).numNodes == 0:
                
                #Mesh a GS that is partially outside the RVE
                Remesh_partial_GS(modelName, gs_part_str, selectedElementCode, eeMeshSize)
                
    #Generate the mesh for the matrix
    Generate_Matrix_Mesh(modelName, matrixName, selectedElementCode)
        
    #This seems to be required by Abaqus
    mdb.models[modelName].rootAssembly.regenerate()


#Mesh GS
def Mesh_GS(model, code, meshSize, gs_part_str):
    
    GS_model = mdb.models[model]
    
    #ElemType(elemCode=code[0], elemLibrary=STANDARD, secondOrderAccuracy=OFF, distortionControl=DEFAULT), 
    #kinematicSplit=AVERAGE_STRAIN, hourglassControl=STIFFNESS, distortionControl=ON, lengthRatio=0.1
    GS_model.parts[gs_part_str].setElementType(
        elemTypes=(
            ElemType(elemCode=code[0], elemLibrary=elementLibrary[selectedLibrary], secondOrderAccuracy=OFF, distortionControl=DEFAULT), 
            ElemType(elemCode=code[1], elemLibrary=elementLibrary[selectedLibrary]), 
            ElemType(elemCode=code[2], elemLibrary=elementLibrary[selectedLibrary])), 
        regions=(GS_model.parts[gs_part_str].cells, ))

    GS_model.parts[gs_part_str].seedPart(deviationFactor=0.1, minSizeFactor=0.1, size=meshSize)
    GS_model.parts[gs_part_str].generateMesh()
       
#Mesh a GS that is partially outside the RVE
def Remesh_partial_GS(modelName, gs_part_str, selectedElementCode, eeMeshSize):
                
    #Set tetrahedron elements in mesh controls
    mdb.models[modelName].parts[gs_part_str].setMeshControls(
        elemShape=TET, regions=mdb.models[modelName].parts[gs_part_str].cells.getSequenceFromMask(('[#1 ]',), ), technique=FREE)

    #Set the element type
    mdb.models[modelName].parts[gs_part_str].setElementType(
        elemTypes=(
            ElemType(elemCode=selectedElementCode[0], elemLibrary=elementLibrary[selectedLibrary], secondOrderAccuracy=OFF, distortionControl=DEFAULT),
            ElemType(elemCode=selectedElementCode[1], elemLibrary=elementLibrary[selectedLibrary], secondOrderAccuracy=OFF, distortionControl=DEFAULT),
            ElemType(elemCode=selectedElementCode[2], elemLibrary=elementLibrary[selectedLibrary], secondOrderAccuracy=OFF, distortionControl=DEFAULT)),
        regions=(mdb.models[modelName].parts[gs_part_str].cells.getSequenceFromMask(('[#1 ]',), ), ))

    #Seed the part
    mdb.models[modelName].parts[gs_part_str].seedPart(deviationFactor=0.1, minSizeFactor=0.1, size=eeMeshSize)

    #Generate mesh
    mdb.models[modelName].parts[gs_part_str].generateMesh()

#Mesh the matrix
def Generate_Matrix_Mesh(modelName, matrixName, selectedElementCode):
                
    #Generate the mesh for the matrix
    mdb.models[modelName].parts[matrixName].setElementType(
        elemTypes=(
            ElemType(elemCode=selectedElementCode[0], elemLibrary=elementLibrary[selectedLibrary], secondOrderAccuracy=OFF, distortionControl=DEFAULT),
            ElemType(elemCode=selectedElementCode[1], elemLibrary=elementLibrary[selectedLibrary]),
            ElemType(elemCode=selectedElementCode[2], elemLibrary=elementLibrary[selectedLibrary])),
        regions=(mdb.models[modelName].parts[matrixName].cells.getSequenceFromMask(('[#1 ]', ), ), ))
    mdb.models[modelName].parts[matrixName].seedPart(deviationFactor=0.1, minSizeFactor=0.1, size=matrixMeshSize)
    mdb.models[modelName].parts[matrixName].generateMesh()
        
        
#Sets needed to create the PBC equations
def Create_Set_for_PBC(model, matrixName, P0, Lxyz, corner):

    modelRoot = mdb.models[model].rootAssembly
    
    #Calculate half the lengths of the RVE along each direction
    half_x = Lxyz[0]*0.5
    half_y = Lxyz[1]*0.5
    half_z = Lxyz[2]*0.5
    
    #Set containing the 6 faces of the RVE (temperature will be applied to this set)
    modelRoot.Set(
        faces=modelRoot.instances[matrixName + '-1'].faces.findAt(
            ((corner[0], half_y, half_z),),
            ((P0[0], half_y, half_z),),
            ((half_x, corner[1], half_z),),
            ((half_x, P0[1], half_y),),
            ((half_x, half_y, corner[2]),),
            ((half_x, half_y, P0[2]),),),
        name='External_Faces')

    #If RPs are needed, then this sould be enabled
    #Creating reference points further than the RVE length and taking its ID
    RP_X_id = modelRoot.ReferencePoint(point=(corner[0]+Lxyz[0], half_y, half_z)).id
    RP_Y_id = modelRoot.ReferencePoint(point=(half_x, corner[1]+Lxyz[1], half_z)).id
    RP_Z_id = modelRoot.ReferencePoint(point=(half_x, half_y, corner[2]+Lxyz[2])).id

    #Creating a set for each reference point (used in the PBC equations)
    modelRoot.Set(name='RP_X', referencePoints=(modelRoot.referencePoints[RP_X_id],))
    modelRoot.Set(name='RP_Y', referencePoints=(modelRoot.referencePoints[RP_Y_id],))
    modelRoot.Set(name='RP_Z', referencePoints=(modelRoot.referencePoints[RP_Z_id],))

    #Set for each face of the RVE (used in the PBC equations)
    modelRoot.Set(faces=
        modelRoot.instances[matrixName + '-1'].faces.findAt(((corner[0], half_y, half_z), )), name='X_positive')
    modelRoot.Set(faces=
        modelRoot.instances[matrixName + '-1'].faces.findAt(((P0[0], half_y, half_z), )), name='X_negative')
    modelRoot.Set(faces=
        modelRoot.instances[matrixName + '-1'].faces.findAt(((half_x, corner[1], half_z), )), name='Y_positive')
    modelRoot.Set(faces=
        modelRoot.instances[matrixName + '-1'].faces.findAt(((half_x, P0[1], half_z), )), name='Y_negative')
    modelRoot.Set(faces=
        modelRoot.instances[matrixName + '-1'].faces.findAt(((half_x, half_y, corner[2]), )), name='Z_positive')
    modelRoot.Set(faces=
        modelRoot.instances[matrixName + '-1'].faces.findAt(((half_x, half_y, P0[2]), )), name='Z_negative')

#Equations of the periodic boundary conditions
def PBC_Equations(model):

    modelRoot = mdb.models[model]

    #Equations with the same reference point (RP_#) terms will be added
    #E.g.: 1*X_positive + 0.5*RP_X + (-1*X_negative) + 0.5*RP_X = X_Positive - X_Negative + RP_X = 0
    #If RPs are needed, then this sould be enabled

    modelRoot.Equation(name='EQ_X_Positive', terms=((1.0, 'X_positive', 1), (-0.5, 'RP_X', 1)))
    modelRoot.Equation(name='EQ_X_Negative', terms=((-1.0, 'X_negative', 1), (-0.5, 'RP_X', 1)))

    modelRoot.Equation(name='EQ_Y_Positive', terms=((1.0, 'Y_positive', 2), (-0.5, 'RP_Y', 2)))
    modelRoot.Equation(name='EQ_Y_Negative', terms=((-1.0, 'Y_negative', 2), (-0.5, 'RP_Y', 2)))

    modelRoot.Equation(name='EQ_Z_Positive', terms=((1.0, 'Z_positive', 3), (-0.5, 'RP_Z', 3)))
    modelRoot.Equation(name='EQ_Z_Negative', terms=((-1.0, 'Z_negative', 3), (-0.5, 'RP_Z', 3)))

#GS vertices for instances (after translations and rotation)
def GS_Instances_NodeSet(model, gs_part_str, gs_inst_str):

    GS_assembly = mdb.models[model].rootAssembly

    #Get all vertices from a graphene sheet
    GS_vertex = GS_assembly.instances[gs_inst_str].vertices
        
    GS_vertices = GS_vertex.getClosest(coordinates=(
        (GS_assembly.getCoordinates(GS_vertex[0])),
        (GS_assembly.getCoordinates(GS_vertex[1])),
        (GS_assembly.getCoordinates(GS_vertex[2])),
        (GS_assembly.getCoordinates(GS_vertex[3])),
        (GS_assembly.getCoordinates(GS_vertex[4])),
        (GS_assembly.getCoordinates(GS_vertex[5])),
        (GS_assembly.getCoordinates(GS_vertex[6])),
        (GS_assembly.getCoordinates(GS_vertex[7])),)) 

    GS_assembly.Set(name=gs_part_str+'_N-0', vertices=(GS_vertex.findAt((((GS_vertices[4][1])),)),))
    GS_assembly.Set(name=gs_part_str+'_N-1', vertices=(GS_vertex.findAt((((GS_vertices[1][1])),)),))
    GS_assembly.Set(name=gs_part_str+'_N-2', vertices=(GS_vertex.findAt((((GS_vertices[0][1])),)),))
    GS_assembly.Set(name=gs_part_str+'_N-3', vertices=(GS_vertex.findAt((((GS_vertices[6][1])),)),))
    GS_assembly.Set(name=gs_part_str+'_N-4', vertices=(GS_vertex.findAt((((GS_vertices[5][1])),)),))
    GS_assembly.Set(name=gs_part_str+'_N-5', vertices=(GS_vertex.findAt((((GS_vertices[2][1])),)),))
    GS_assembly.Set(name=gs_part_str+'_N-6', vertices=(GS_vertex.findAt((((GS_vertices[3][1])),)),))
    GS_assembly.Set(name=gs_part_str+'_N-7', vertices=(GS_vertex.findAt((((GS_vertices[7][1])),)),))

#This function creates two sets for two of the vertices of the matrix (sample), where each set has only one node
#These nodes are on the diagonal of the cuboid that defines the matrix (sample)
def Create_Sets_for_Matrix(model, P0, corner, matrixName):
    
    #Create the set of the lower left corner
    mdb.models[model].rootAssembly.Set(
        vertices=mdb.models[model].rootAssembly.instances[matrixName+'-1'].vertices.getByBoundingSphere(center=P0, radius=0.001),
        name='Matrix0')
   
    #Create the set of the opposite corner
    mdb.models[model].rootAssembly.Set(
        vertices=mdb.models[model].rootAssembly.instances[matrixName+'-1'].vertices.getByBoundingSphere(center=corner, radius=0.001),
        name='Matrix1')

def Data_Check_And_Remesh_GSs(modelName, jobName, selectedElementCode, eeMeshSize):
    
    #Perform data check
    start = time.time()
    mdb.jobs[jobName].submit(consistencyChecking=OFF, datacheckJob=True)
    mdb.jobs[jobName].waitForCompletion()
    
    #Name of lock file
    lckFile = jobName+'.lck'
    
    #Wait for lck file to dissapear
    while os.path.exists(lckFile):
        #print('waiting...')
        time.sleep(5)
        
    plog("Data check on input file: {} secs.\n".format(time.time()-start))
    
    #Check if there is a set for elements with zero volume
    try: 
        
        #Open database
        odb = openOdb(path=jobName+'.odb')
        #print(odb.rootAssembly.elementSets['ErrElemVolSmallNegZero'])
        #print(odb.rootAssembly.elementSets['ErrElemVolSmallNegZero'].instanceNames)
        
        #Get the instances with zero elements if any
        instanceNames = odb.rootAssembly.elementSets['ErrElemVolSmallNegZero'].instanceNames
            
        #Close odb file
        odb.close()
        
        #Iterate over the instances with zero volume elements 
        for inst in instanceNames:
            #print(inst)
            #print(inst[0:len(inst)-2])
            
            #Get part name
            gs_part_str = inst[0:len(inst)-2]
            plog('Remeshing part '+gs_part_str+' due to zero volume elements\n')
            
            #Mesh a GS that is partially outside the RVE
            Remesh_partial_GS(modelName, gs_part_str, selectedElementCode, eeMeshSize)
            
        #This seems to be required by Abaqus to update mesh
        mdb.models[modelName].rootAssembly.regenerate()
        
        #Try the data check on the Abaqus input file again 
        start = time.time()
        mdb.jobs[jobName].submit(consistencyChecking=OFF, datacheckJob=True)
        mdb.jobs[jobName].waitForCompletion()
        
        #Wait for lck file to dissapear
        while os.path.exists(lckFile):
            #print('waiting...')
            time.sleep(5)
            
        plog("Second data check on input file: {} secs.\n".format(time.time()-start))
        
    except:
        plog('Set ErrElemVolSmallNegZero not found.\n')


##############################################################################################
##############################################################################################
##############################################################################################
#####################################---MAIN PROGRAM---#######################################


######################################---OPENING CSV---#######################################

start0 = time.time()

#with open(os.path.join(os.path.dirname('__file__'), csv_gnpFile + '.csv')) as file_gnp:
with open(csv_gnpFile) as file_gnp:
    data_gnp = list(csv.reader(file_gnp, quoting=csv.QUOTE_NONNUMERIC))

#with open(os.path.join(os.path.dirname('__file__'), csv_geomFile + '.csv')) as file_geom:
with open(csv_geomFile) as file_geom:
    #data_geom = list(csv.reader(file_geom, quoting=csv.QUOTE_NONNUMERIC))
    (P0, Lxyz) = list(csv.reader(file_geom, quoting=csv.QUOTE_NONNUMERIC))

#Get the corner of the RVE opposite to P0
corner = (P0[0]+Lxyz[0], P0[1]+Lxyz[1], P0[2]+Lxyz[2])

#Calculate the number of GNPs
N_GSs = len(data_gnp)

#Name of the job to be used based on its parameters
#GS-'Number of GSs in the RVE'
#EPS-'Number of elements per side'
#MR-'Mesh ratio' (EmbeddedMesh/HostMesh)x100
jobName = 'GS-'+str(N_GSs)+'_EPS-'+str(elementsPerSide)+'_MR-'+str(int(100*meshRatio))

#Name of the file to save print messages
print_file = jobName + '.txt'
pfile = open(print_file, "a")
plog('####################  START  ####################\n')
plog('There are ' + str(N_GSs) + ' GSs inside the RVE.\n')
plog('Reading csv files time: {}\n'.format(time.time()-start0))

#Number of GS laying inside/outside the RVE
indexOutside = []
indexInside = []

#Select an element type based on the choice of geometric order
selectedElementCode = Select_Elemet_Type(selectedGeomOrder)

#Mesh size for the matrix (um)
matrixMeshSize = Lxyz[1]/elementsPerSide

#Embedded element mesh size
eeMeshSize = matrixMeshSize*meshRatio

####################################---MODEL CREATION---######################################

start = time.time()

#Creating the matrix
Create_Matrix(modelName, sheetSize, matrixName, P0, Lxyz)

#Creating a bounding box
Create_BiggerBox(modelName, sheetSize, biggerBoxName, P0, Lxyz, margin)

#Creating each material
Create_Material(modelName, matrixMaterial, matrixDensity, matrixModulus, matrixPoissonR, matrixExpCoeff, matrixThermConductivity, matrixSpecHeat, matrixElecConductivity)
Create_Material(modelName, fillerMaterial, fillerDensity, fillerModulus, fillerPoissonR, fillerExpCoeff, fillerThermConductivity, fillerSpecHeat, fillerElecConductivity)

#Creating a section for each material
Create_Section(modelName, matrixMaterial)
Create_Section(modelName, fillerMaterial)

#-----------------------------------START: CREATE ASSEMBLY-----------------------------------#

#This line seems to be required by Abaqus to set a Cartesian system
mdb.models[modelName].rootAssembly.DatumCsysByDefault(CARTESIAN)

#Create instance of the matrix
Create_Matrix_Instance(modelName, matrixName, P0)

#Create instance of the bigger box
mdb.models[modelName].rootAssembly.Instance(dependent=ON, name=biggerBoxName + '-1', part=mdb.models[modelName].parts[biggerBoxName])

#Create the box that will be used to cut all GS that are partially outside the RVE
Create_CuttingBox(modelName, matrixName, biggerBoxName, P0, margin)

#Assign section to matrix part
Assign_Section(modelName, matrixMaterial, matrixName)
#mdb.models[modelName].parts[matrixName].SectionAssignment(
#    offset=0.0, offsetField='', offsetType=MIDDLE_SURFACE,
#    region=Region(cells=mdb.models[modelName].parts[matrixName].cells),
#    sectionName=matrixMaterial, thicknessAssignment=FROM_SECTION)

#Create parts and instances for GSs
Create_All_GSs(modelName, fillerMaterial, sheetSize, N_GSs, P0, corner)

end = time.time()
plog("Time for part and instance generation: {}\n".format(end-start))

#------------------------------------END: CREATE ASSEMBLY------------------------------------#

#-----------------------------START: PERIODIC BOUNDARY CONDITIONS----------------------------#

#Defining periodic boundary conditions

#Create sets for each face of the RVE
Create_Set_for_PBC(modelName, matrixName, P0, Lxyz, corner)

#Create the sets for
Create_Sets_for_Matrix(modelName, P0, corner, matrixName)

#Add the equations that define the PBC
PBC_Equations(modelName)

#Add temperature boundary conditions
Add_Boundary_Conditions(modelName, tempApplied)

#print('Periodic boundary conditions have been applied.')

#Set the output request
mdb.models[modelName].fieldOutputRequests['F-Output-1'].setValues(variables=('S', 'E', 'U'))

#------------------------------END: PERIODIC BOUNDARY CONDITIONS-----------------------------#

#---------------------------------------START: MESHING---------------------------------------#
start = time.time()

Generate_Meshes(modelName, matrixName, selectedElementCode, eeMeshSize, N_GSs, indexInside, indexOutside)

end = time.time()
plog("Time for meshing: {}\n".format(end-start))

#Check if files for re-meshing the model are needed
if reMeshModel == 1:
    #Create files that contain the index of each GS depending on its state (inside or outside the RVE)
    textfile_IO = open('GS-'+str(N_GSs)+'_indexOutside.txt', 'w')
    for element in indexOutside:
        textfile_IO.write(str(element)+',')
    textfile_IO.close()

    textfile_II = open('GS-'+str(N_GSs)+'_indexInside.txt', 'w')
    for element in indexInside:
        textfile_II.write(str(element)+',')
    textfile_II.close()
#---------------------------------------END: MESHING----------------------------------------#

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
        #Also remesh GSs with zero volume elements if any
        Data_Check_And_Remesh_GSs(modelName, jobName, selectedElementCode, eeMeshSize)

    if submitJob == 1:
        
        #Submit job and save time of submission
        start = time.time()
        mdb.jobs[jobName].submit(consistencyChecking=OFF)
        
        #Make the Python script to wait for the Abaqus job to finish
        #In this way the script can measure the execution time of the Abaqus simulation
        mdb.jobs[jobName].waitForCompletion()
        
        end = time.time()
        plog("Time for Job execution: {}\n".format(end-start))

end = time.time()        
plog("Time for Abaqus model: {}\n".format(end-start0))
pfile.close()
