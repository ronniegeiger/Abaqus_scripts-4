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
#Name of the .csv file with the GS parameters
csv_gnpFile = 'gnp_data.csv'
#Name of the .csv file with the RVE geometry parameters
csv_geomFile = 'sample_geom.csv'
#Name of the .csv file with the coordinates of CNT points
coord_file = 'cnt_coordinates.csv'
#Name of the .csv file with total number of CNTs, number of points per CNT and radius of each CNT
struct_file = 'cnt_struct.csv'

###################### Mesh flags and variables
#Number of elements per side of the RVE
elementsPerSide = 20
#EmbeddedMesh/HostMesh ratio
meshRatio = 0.75
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
saveInputFile = 1
#Perform data check on input file
#If flag saveInputFile is set to 0, this flag is ignored
# 0 = no
# 1 = yes
dataCheck = 0
#Submit job: 
# 0 = no
# 1 = yes
submitJob = 0

###################### Restart variables
#Restart flag
#0 = no restart
#1 = restart
restartFlag = 1

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
dispX = 0.09
dispY = 0.0
dispZ = 0.0

###################### Matrix properties
#Define the name of the matrix
matrixMaterial = 'Polypropylene'
#Matrix name
matrixName = 'MATRIX'
#RVE name
rveName = 'RVE'
#String for the host set (i.e., the matrix)
strHost = 'HOST_SET'
#Define the mass density (kg/m3)
matrixDensity = 0.9e-9
#Define the elastic modulus (GPa)
matrixModulus = 3.0
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

### GSs
#Define the name of the filler
gsMaterial = 'Graphene'
#Define the mass density (kg/m3)
gsDensity = 2200
#Define the elastic modulus (GPa) - 1000
gsYoungsModulus = 1000
#Define the Poisson ratio - 0.165
gsPoissonR = 0.149
#Define de coefficient of thermal expansion (e-5 1/C)
gsExpCoeff = 0.5
#Define the electrical conductivity (S/m)
gsElecConductivity = 10e7
#Define the thermal conductivity (W/m*K) - Because of Abaqus needs, but doesn't affect the thermo-mechanical simulation (could be = 1)
gsThermConductivity = 3000
#Define the specific heat (J/mol*K) - Because of Abaqus needs, but doesn't affect the thermo-mechanical simulation (could be = 1)
gsSpecHeat = 7
#Maximum dimensions
lxMax = 3.0
lyMax = 3.0
tMax = 0.03
margin = 1.01*sqrt(lxMax*lxMax + lyMax*lyMax + tMax*tMax)

### CNTs
#Define the name and section of the filler
cntMaterial = 'CNT_MAT'
cntSection = 'CNT_SEC'
#Define the mass density (kg/m3)
cntDensity = 2.2e-9
#Define the elastic modulus (GPa)
cntYoungsModulus = 1000
#Define the Poisson ratio
cntPoissonR = 0.165
#Define de coefficient of thermal expansion (e-5 1/C)
cntExpCoeff = 0.5
#Define the thermal conductivity (W/m*K)
cntThermConductivity = 3000
#Define the specific heat (J/mol*K)
cntSpecHeat = 7
#Define the electrical conductivity (S/m)
cntElecConductivity = 1e7
#Maximum CNT radius
cnt_rad_max = 0.015

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

#Cosine of 45 degrees (PI/4)
cos45 = 0.7071067812

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

#This function generates a string for the node set of CNT i
def cnt_string_node_set(cnt_i):
    #Return the string in the format CNT-nodes-cnt_i
    return 'CNT-%d-NODES' %(cnt_i)

def ee_string(filler, i):
    return 'EE-%s-%d'%(filler, i)

######################################---MATH FUNCTIONS---########################################

#This function generates a rotation matrix
def rotation_matrix(P1, P2):

    #Get the components of the vector in the direction of P1P2
    ux = P2[0] - P1[0]
    uy = P2[1] - P1[1]
    uz = P2[2] - P1[2]

    #Squared quantities
    ux2 = ux*ux
    uy2 = uy*uy

    #Calculate the length of the vector u
    u_length = sqrt(ux2 + uy2 + uz*uz);

    #This quantity is used three times:
    quantity = sqrt(ux2 + uy2);

    #Calculate the trigonometric functions of the angles theta and phi
    cos_phi = ux/quantity;
    sin_phi = uy/quantity;
    cos_theta = uz/u_length;
    sin_theta = quantity/u_length;

    #Generate the components of the rotation matrix
    R11 = cos_phi*cos_theta
    R12 = -sin_phi
    R13 = cos_phi*sin_theta

    R21 = sin_phi*cos_theta
    R22 = cos_phi
    R23 = sin_phi*sin_theta

    R31 = -sin_theta
    R33 = cos_theta

    #Create a 'matrix' using the components found above
    R = ((R11, R12, R13),(R21, R22, R23),(R31, 0.0, R33))

    #Return the rotation matrix
    return R

#This function retunrs a unit vector going from P1 towards P2
def get_unit_vector(P1, P2):

    #Get the vector going from P1 towards P2
    v = (P2[0]-P1[0], P2[1]-P1[1], P2[2]-P1[2])

    #Calculate the length of vector v
    length = sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2])

    #Return v as a unit vector
    return (v[0]/length, v[1]/length, v[2]/length)

#This function makes the input vector a unit vector
def make_unit(V):

    #Calculate length of V
    length = sqrt(V[0]*V[0] + V[1]*V[1] + V[2]*V[2])

    #Return V as a unit vector
    return (V[0]/length, V[1]/length, V[2]/length)

#Cross product of two vectors
def cross(v1, v2):

    #Calculate the compoenents of the cross product
    x = v1[1]*v2[2] - v1[2]*v2[1]
    y = -v1[0]*v2[2] + v1[2]*v2[0]
    z = v1[0]*v2[1] - v1[1]*v2[0]

    #Return the cross product
    return (x,y,z)

#Dot product of two vectors
def dot(v1, v2):
    return (v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2])

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
def Create_RVE(model, sheetsz, rveName, P0, Lxyz):

    #Create the matrix geometry
    #The matrix is the extended RVE
    mdb.models[model].ConstrainedSketch(name='__profile__', sheetSize=sheetsz)
    #Define the two corners of a rectangle on the xy plane
    mdb.models[model].sketches['__profile__'].rectangle(
        point1=(P0[0], P0[1]), point2=(P0[0]+Lxyz[0], P0[1]+Lxyz[1]))
    #Name the part
    mdb.models[model].Part(dimensionality=THREE_D, name=rveName, type=DEFORMABLE_BODY)
    
    #Use the length along the z-direction as the extrusion depth
    mdb.models[model].parts[rveName].BaseSolidExtrude(
        depth=Lxyz[2], sketch=mdb.models[model].sketches['__profile__'])
    
    #Delete the sketch
    del mdb.models[model].sketches['__profile__']
    
def Create_Matrix(modelName, P0, Lxyz, Lxyz_ext, matrixName, matrixMeshSize, halfMatrixMeshSize):
    
    #Create a sketch
    mdb.models[modelName].ConstrainedSketch(name='__profile__', sheetSize=200.0)
    
    #Define the two corners of a rectangle on the xy plane
    mdb.models[modelName].sketches['__profile__'].rectangle(
        point1=(P0[0], P0[1]), 
        point2=(P0[0]+Lxyz_ext[0], P0[1]+Lxyz_ext[1]))
        
    #Name the part
    mdb.models[modelName].Part(
        dimensionality=THREE_D, 
        name=matrixName, 
        type=DEFORMABLE_BODY)
    
    #Use the length along the z-direction as the extrusion depth
    mdb.models[modelName].parts[matrixName].BaseSolidExtrude(
        depth=Lxyz_ext[2], 
        sketch=mdb.models[modelName].sketches['__profile__'])
        
    #Delete the sketch
    del mdb.models[modelName].sketches['__profile__']
    
    #Partition cell by point normal

    #First cut
    #Reference point for cutting the cell
    P_cut = [P0[0]+Lxyz[0]+matrixMeshSize, P0[1], Lxyz_ext[2] ]
    #Get reference edge for first and second cut
    normal_edge = mdb.models[modelName].parts[matrixName].edges.findAt(P_cut, )
    #Cut all cells
    mdb.models[modelName].parts[matrixName].PartitionCellByPlanePointNormal(
        cells=mdb.models[modelName].parts[matrixName].cells, 
        normal=normal_edge, 
        point=P_cut)
    
    #Second cut
    #Update reference point for cutting the cell
    P_cut[0] = P0[0]+matrixMeshSize
    #Cut cell found using reference point
    mdb.models[modelName].parts[matrixName].PartitionCellByPlanePointNormal(
        cells=mdb.models[modelName].parts[matrixName].cells.findAt(P_cut, ), 
        normal=normal_edge, 
        point=P_cut)
    
    #Third cut
    #Update reference point for cutting the cell
    P_cut = [P0[0], P0[1]+matrixMeshSize, Lxyz_ext[2] ]
    #Datum point for debugging
    #mdb.models[modelName].parts[matrixName].DatumPointByCoordinate(P_cut)
    #Update reference edge
    normal_edge = mdb.models[modelName].parts[matrixName].edges.findAt(P_cut, )
    #print(normal_edge)
    #Cut all cells
    mdb.models[modelName].parts[matrixName].PartitionCellByPlanePointNormal(
        cells=mdb.models[modelName].parts[matrixName].cells, 
        normal=normal_edge, 
        point=P_cut)
    
    #Fourth cut
    #Update reference point for cutting the cell
    P_cut[1] = P0[1]+Lxyz[1]+matrixMeshSize
    #Datum point for debugging
    #mdb.models[modelName].parts[matrixName].DatumPointByCoordinate(P_cut)
    #Find cells to be cut
    cells_tmp = mdb.models[modelName].parts[matrixName].cells.getByBoundingBox(
        P0[0]-halfMatrixMeshSize, P0[1]+halfMatrixMeshSize, -halfMatrixMeshSize, 
        P0[0]+Lxyz_ext[0]+halfMatrixMeshSize, P0[1]+Lxyz_ext[1]+halfMatrixMeshSize, Lxyz_ext[2]+halfMatrixMeshSize)
    #print('cells_tmp ',cells_tmp, ' len=', len(cells_tmp))
    #print(cells_tmp[0])
    #print(cells_tmp[1])
    #Update reference edge
    normal_edge = mdb.models[modelName].parts[matrixName].edges.findAt(P_cut, )
    #print(normal_edge)
    #Cut cells
    mdb.models[modelName].parts[matrixName].PartitionCellByPlanePointNormal(
        cells=(cells_tmp[0], cells_tmp[1], cells_tmp[2]),
        normal=normal_edge, 
        point=P_cut)
    
    #Fifth cut
    #Update reference point for cutting the cell
    P_cut = [P0[0], P0[1], Lxyz[2]+matrixMeshSize ]
    #Update reference edge
    normal_edge = mdb.models[modelName].parts[matrixName].edges.findAt(P_cut, )
    #Cut all cells
    mdb.models[modelName].parts[matrixName].PartitionCellByPlanePointNormal(
        cells=mdb.models[modelName].parts[matrixName].cells, 
        normal=normal_edge, 
        point=P_cut)
    
    #Sixth cut
    #Update reference point for cutting the cell
    P_cut[2] = matrixMeshSize
    #Update reference edge
    normal_edge = mdb.models[modelName].parts[matrixName].edges.findAt(P_cut, )
    #Find cells to be cut
    cells_tmp = mdb.models[modelName].parts[matrixName].cells.getByBoundingBox(
        P0[0]-halfMatrixMeshSize, P0[1]-halfMatrixMeshSize, -halfMatrixMeshSize, 
        P0[0]+Lxyz_ext[0]+halfMatrixMeshSize, P0[1]+Lxyz_ext[1]+halfMatrixMeshSize, Lxyz_ext[2]-halfMatrixMeshSize)
    #Cut cells using a bounding box
    mdb.models[modelName].parts[matrixName].PartitionCellByPlanePointNormal(
        cells=tuple(cells_tmp), 
        normal=normal_edge, 
        point=P_cut)

#Create a box bigger than the RVE, which is used to trim the GNPs
def Create_BiggerBox(model, sheetsz, part, P0, Lxyz, margin):
    
    #The length of the bigger box along each direction is the same as the RVE plus the margin
    lengthBiggerBox_x = Lxyz[0]+2.0*margin
    lengthBiggerBox_y = Lxyz[1]+2.0*margin
    lengthBiggerBox_z = Lxyz[2]+2.0*margin

    #Create the geometry for the bigger box
    mdb.models[model].ConstrainedSketch(name='__profile__', sheetSize=sheetsz)
    #Note that coordinates of the min point will be (P0[0]-margin, P0[1]-margin, 0)
    mdb.models[model].sketches['__profile__'].rectangle(
        point1=(P0[0]-margin, P0[1]-margin), 
        point2=(P0[0]-margin+lengthBiggerBox_x, P0[1]-margin+lengthBiggerBox_y))
    mdb.models[model].Part(dimensionality=THREE_D, name=part, type=DEFORMABLE_BODY)

    mdb.models[model].parts[part].BaseSolidExtrude(
        depth=lengthBiggerBox_z, sketch=mdb.models[model].sketches['__profile__'])
    del mdb.models[model].sketches['__profile__']

def Create_Matrix_Instance(modelName, matrixName, P0, matrixMeshSize):

    #Instance name of matrix
    str_mat_inst = matrixName + '-1'

    #Create an instance of the matrix
    mdb.models[modelName].rootAssembly.Instance(dependent=ON, name=str_mat_inst, part=mdb.models[modelName].parts[matrixName])

    #Translate instance
    #i.e: endpoint - starting point = (P0[0], P0[1], P0[2]) - (P0[0]+matrixMeshSize, P0[1]+matrixMeshSize, matrixMeshSize)
    mdb.models[modelName].rootAssembly.translate(
        instanceList=(str_mat_inst, ),
        vector=(-matrixMeshSize, -matrixMeshSize, P0[2]-matrixMeshSize))

def Create_RVE_Instance(modelName, rveName, P0):

    #Instance name of RVE
    rveIntance = rveName + '-1'

    #Create an instance of the RVE
    mdb.models[modelName].rootAssembly.Instance(
        dependent=ON, 
        name=rveIntance, 
        part=mdb.models[modelName].parts[rveName])

    #Check if the instance needs to bre translated along the z-axis
    if abs(P0[2]) > Zero:

        #z-coordinate of "min point" of RVE is non-zero
        #Current coordinates of min point are (P0[0], P0[1], 0)
        #Move that corner to (P0[0], P0[1], P0[2])
        #i.e: endpoint - starting point = (P0[0], P0[1], P0[2]) - (P0[0], P0[1], 0)
         mdb.models[modelName].rootAssembly.translate(
             instanceList=(rveIntance, ),
             vector=(0.0, 0.0, P0[2]))

#Create the hollow box that is used to cut all GS that are partially outside the RVE
def Create_CuttingBox(model, rveName, partBox, P0, margin):
    
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
    #Both BiggerBox and RVE are not needed, so they are suppressed after creating the Cutter
    mdb.models[model].rootAssembly.InstanceFromBooleanCut(
        cuttingInstances=(mdb.models[model].rootAssembly.instances[rveName + '-1'], ),
        instanceToBeCut=mdb.models[model].rootAssembly.instances[partBox + '-1'],
        name='CUTTER',
        originalInstances=SUPPRESS)

#This function generates all CNT parts
def CNT_Parts_All(model, N_CNTs, cnt_struct, cnt_coords):
    
    #Number of accumulated points
    acc_pts = 0
    
    #Variables to check progress
    check_step = 0.1
    frac_thres = 0.1
    
    #Get the starting time
    start_cnts = time.time()

    #Iterate over the number of CNTs
    for cnt_i in range(1, N_CNTs+1):

        #Number of points in CNTi and its radius
        N_p = int(cnt_struct[cnt_i][0])
        rad = cnt_struct[cnt_i][1]

        #Create all CNTs, one part per CNT
        #The first point is given by the accumulated number of points
        #The last point is the accumulated number of points plus the number of points of CNTi minus 1
        CNT_Part(model, cnt_i, rad, acc_pts, acc_pts+N_p-1, cnt_coords)
        
        #Calculate the fraction of generated CNT parts
        frac = float(cnt_i)/float(N_CNTs)
        
        #Check if more than the fraction threshold has been generated
        if frac >= frac_thres:
            
            #Update the threshold
            while frac >= frac_thres:
                
                #Increase the threshold by a step
                frac_thres += check_step
                
            #Send message with percentage generated
            plog('   CNT parts generated: {} % ({} secs)\n'.format( int( (frac_thres - check_step)*100 ), time.time()-start_cnts))

        #Increase the number of accumulated points
        acc_pts += N_p

#This function creates a CNT in Part module
def CNT_Part(model, cnt_i, cnt_rad, cnt_start, cnt_end, cnt_coords):
	
    #Get the string for the CNT
    str_part = string_part('CNT', cnt_i)

    #Create a point to be able to generate the edges that will make the CNT
    mdb.models[model].Part(dimensionality=THREE_D, name=str_part, type=DEFORMABLE_BODY)
    mdb.models[model].parts[str_part].ReferencePoint(point=(0.0, 0.0, 0.0))

    #print("cnt_start=",cnt_start," cnt_end=",cnt_end)

    #Create edges, iterate over all points of the CNT
    Generate_Edges(model, cnt_start, cnt_end, cnt_coords, str_part)

    #Sweep an octagon along the edges
    Generate_Sweep(model, cnt_i, cnt_coords[cnt_end-1], cnt_coords[cnt_end], cnt_rad, cnt_start, cnt_end, str_part)

    #Delete the initial point as it is not used anymore
    del mdb.models[model].parts[str_part].features['RP']

#This function generates all edges of a CNT
def Generate_Edges(model, cnt_start, cnt_end, cnt_coords, str_part):
    
    #Create edges, iterate over all points of the CNT
    for i in range(cnt_start, cnt_end):
    #for i in range(cnt_start, cnt_start+11):

        #First point of the CNT segment
        P1 = cnt_coords[i]
        #Second point of the CNT segment
        P2 = cnt_coords[i+1]

        #Generate the segment
        mdb.models[model].parts[str_part].WirePolyLine(mergeType=IMPRINT, meshable= ON, points=((P1, P2), ))

        #Generate the name of the wire
        #str_wire = 'Wire-%d-Set-1' %(i-cnt_start+1)

        #Name the wire
        #mdb.models[model].parts[str_part].Set(edges=
        #	mdb.models[model].parts[str_part].edges.getSequenceFromMask(('[#1 ]', ), ), name=str_wire)

#This function generates the sweeped part
def Generate_Sweep(model, cnt_i, P1, P2, cnt_rad, cnt_start, cnt_end, str_part):

    #Profile name is used several times, so create a variable to modify it easyly
    #in case it is needed
    profile_str = '__profile__'
    
    #Generate the octagon that will be swept through the CNT segments
    Generate_Sweep_Profile(model, P1, P2, cnt_rad, str_part, profile_str)
    
    #Perform Sweep
    try:
        Sweep_Command_profileNormal_OFF(model, str_part, profile_str, cnt_start, cnt_end)
    except:
        
        plog('\nCould not generate sweep of part '+str_part+'\n')
        plog('cnt_end= {} cnt_start={}\n'.format(cnt_end, cnt_start))
        plog('Trying again with profileNormal=ON\n\n')
        Sweep_Command_profileNormal_ON(model, str_part, profile_str, cnt_start, cnt_end)

    #Delete sketch
    del mdb.models[model].sketches[profile_str]

#This function generates the octagon that will be swet through the CNT segmetns
def Generate_Sweep_Profile(model, P1, P2, cnt_rad, str_part, profile_str):
    #The sketching plane is perpendicular to the last segment
    #The last segment has the points P1 and P2, so the z-axis of the sketch plane is 
    #aligned to this last segment and goes in the direction from P1 to P2
    #A transformation matrix is generated to align the z-axis to the last segment
    R = rotation_matrix(P1, P2)

    #Create the sketching plane using the rotation matrix and the last point in the CNT
    mdb.models[model].ConstrainedSketch(gridSpacing=0.001, name=profile_str, sheetSize=0.076, transform=(
        R[0][0], R[0][1], R[0][2],
        R[1][0], R[1][1], R[1][2], 
        R[2][0], 0.0, R[2][2], 
        P2[0], P2[1], P2[2]))
    mdb.models[model].sketches[profile_str].sketchOptions.setValues( decimalPlaces=5)
    mdb.models[model].sketches[profile_str].ConstructionLine(point1=(-0.038, 0.0), point2=(0.038, 0.0))
    mdb.models[model].sketches[profile_str].ConstructionLine(point1=(0.0,-0.038), point2=(0.0, 0.038))
    mdb.models[model].parts[str_part].projectReferencesOntoSketch(
        filter=COPLANAR_EDGES, sketch=mdb.models[model].sketches[profile_str])

    #Calculate the radius multiplied by cos(PI/4)
    new_rad = cnt_rad*cos45

    #Construct a regular octagon
    #Vertex 1-2
    mdb.models[model].sketches[profile_str].Line(point1=(0.0, cnt_rad), point2=(new_rad, new_rad))
    #Vertex 2-3
    mdb.models[model].sketches[profile_str].Line(point1=(new_rad, new_rad), point2=(cnt_rad, 0.0))
    #Vertex 3-4
    mdb.models[model].sketches[profile_str].Line(point1=(cnt_rad, 0.0), point2=(new_rad, -new_rad))
    #Vertex 4-5
    mdb.models[model].sketches[profile_str].Line(point1=(new_rad,-new_rad), point2=(0.0, -cnt_rad))
    #Vertex 5-6
    mdb.models[model].sketches[profile_str].Line(point1=(0.0, -cnt_rad), point2=(-new_rad, -new_rad))
    #Vertex 6-7
    mdb.models[model].sketches[profile_str].Line(point1=(-new_rad, -new_rad), point2=(-cnt_rad, 0.0))
    #Vertex 7-8
    mdb.models[model].sketches[profile_str].Line(point1=(-cnt_rad, 0.0), point2=(-new_rad, new_rad))
    #Vertex 8-1
    mdb.models[model].sketches[profile_str].Line(point1=(-new_rad, new_rad), point2=(0.0, cnt_rad))

def Sweep_Command_profileNormal_OFF(model, str_part, profile_str, cnt_start, cnt_end):

    #Select the last edge, which has number equal to (number edges-1)
    #The number of edges is equal to (number of points-1)
    #The number of points is (cnt_end-cnt_start+1)
    #Then, the las edge has number ((cnt_end-cnt_start+1)-1-1)=(cnt_end-cnt_start-1)
    mdb.models[model].parts[str_part].SolidSweep(
        path=mdb.models[model].parts[str_part].edges, 
        profile=mdb.models[model].sketches[profile_str], 
        sketchOrientation=RIGHT, 
        sketchUpEdge=mdb.models[model].parts[str_part].edges[cnt_end-cnt_start-1])

def Sweep_Command_profileNormal_ON(model, str_part, profile_str, cnt_start, cnt_end):

    #Select the last edge, which has number equal to (number edges-1)
    #The number of edges is equal to (number of points-1)
    #The number of points is (cnt_end-cnt_start+1)
    #Then, the las edge has number ((cnt_end-cnt_start+1)-1-1)=(cnt_end-cnt_start-1)
    mdb.models[model].parts[str_part].SolidSweep(
        path=mdb.models[model].parts[str_part].edges, 
        profile=mdb.models[model].sketches[profile_str],
        profileNormal=ON, 
        sketchOrientation=RIGHT, 
        sketchUpEdge=mdb.models[model].parts[str_part].edges[cnt_end-cnt_start-1])

#This function generates all CNT parts
def CNT_Parts_Range(modelName, N_CNT_start, N_CNT_end, cnt_struct, cnt_coords, cntSection):
    
    #Number of accumulated points
    acc_pts = Accumulate_points(N_CNT_start, cnt_struct)
    
    #Variables to check progress
    check_step = 0.1
    frac_thres = 0.1
    
    #Calculate the number of CNTs
    N_CNTs = N_CNT_end - N_CNT_start + 1
    
    #Get the starting time
    start_cnts = time.time()

    #Iterate over the number of CNTs
    for cnt_i in range(N_CNT_start, N_CNT_end+1):

        #Number of points in CNTi and its radius
        N_p = int(cnt_struct[cnt_i][0])
        rad = cnt_struct[cnt_i][1]

        #Create all CNTs, one part per CNT
        #The first point is given by the accumulated number of points
        #The last point is the accumulated number of points plus the number of points of CNTi minus 1
        CNT_Part(modelName, cnt_i, rad, acc_pts, acc_pts+N_p-1, cnt_coords)

        #Get the string for the CNT part
        cnt_str = string_part('CNT', cnt_i)

        #Assign the CNT section to cnt_i
        mdb.models[modelName].parts[cnt_str].SectionAssignment(
            offset=0.0, 
            offsetField='', 
            offsetType=MIDDLE_SURFACE, 
            region=Region(cells=mdb.models[modelName].parts[cnt_str].cells),
            sectionName=cntSection, 
            thicknessAssignment=FROM_SECTION)

        #Generate name for CNT instance
        cnt_inst_str = string_instance('CNT', cnt_i)

        #Generate CNT instance
        mdb.models[modelName].rootAssembly.Instance(
            dependent=ON, 
            name=cnt_inst_str,
            part=mdb.models[modelName].parts[cnt_str])
        
        #Calculate the fraction of generated CNT parts and instances
        frac = float(cnt_i-N_CNT_start+1)/float(N_CNTs)
        
        #Check if more than the fraction threshold has been generated
        if frac >= frac_thres:
            
            #Update the threshold
            while frac >= frac_thres:
                
                #Increase the threshold by a step
                frac_thres += check_step
                
            #Send message with percentage generated
            plog('   CNT parts and instances generated: {} % ({} secs)\n'.format( int( (frac_thres - check_step)*100 ), time.time()-start_cnts))

        #Increase the number of accumulated points
        acc_pts += N_p

def Accumulate_points(N_CNT_start, cnt_struct):
    
    #Check the case of the first CNT being 1
    if N_CNT_start == 1:
        #Nothing to accumulate
        return 0
    
    #Initialize the number of accumulated points
    acc_pts = 0
    
    #Iterate over all CNTs
    #For-loop will stop at CNT with number N_CNT_start-1
    for cnt_i in range(1, N_CNT_start):

        #Number of points in CNTi
        N_p = int(cnt_struct[cnt_i][0])

        #Increase the number of accumulated points
        acc_pts += N_p
    
    #return the number of accumulated points
    return acc_pts

def Accumulate_Points(cnt_start, cnt_end, cnt_struct):
    
    #Check the case of the first CNT being 1
    if cnt_end == 1:
        #Nothing to accumulate
        return 0
    
    #Initialize the number of accumulated points
    acc_pts = 0
    
    #Iterate over all CNTs in the range
    #For-loop will stop at CNT with number cnt_end-1
    for cnt_i in range(cnt_start, cnt_end):

        #Number of points in CNTi
        N_p = int(cnt_struct[cnt_i][0])

        #Increase the number of accumulated points
        acc_pts += N_p
    
    #return the number of accumulated points
    return acc_pts
    
#This function deletes a part and generates a new version
def Delete_Imported_Create_New(modelName, cnts_new, cnt_struct, cnt_coords, cntMaterial, arr_acc_pts):
    
    #Get the starting time
    start_cnts = time.time()
    
    #Variable to store the previous CNT in the for-loop
    cnt_prev = 1
    
    #Initialize the number of accumulated points
    acc_pts = 0

    #Iterate over the number of CNTs
    for cnt_i in cnts_new:
        
        #Accumulate points
        acc_pts += Accumulate_Points(cnt_prev, cnt_i, cnt_struct)
        
        #Append element to vector of accumulated points
        arr_acc_pts.append(acc_pts)

        #Number of points in CNTi and its radius
        N_p = int(cnt_struct[cnt_i][0])
        rad = cnt_struct[cnt_i][1]
        
        #Generate part name
        str_part = string_part('CNT', cnt_i)

        #Generate name for CNT instance
        str_inst = string_instance('CNT', cnt_i)
        
        #Delete old sets for cnt_i
        del mdb.models[modelName].rootAssembly.sets[str_part+'-NODES']
        del mdb.models[modelName].rootAssembly.sets['EE-'+str_part]
        
        #Delete old instance
        del mdb.models[modelName].rootAssembly.features[str_inst]
        
        #Delete old part
        del mdb.models[modelName].parts[str_part]
        
        plog('Part '+str_part+' and instance deleted')
        
        #Delete old embedded element constraint
        del mdb.models[modelName].constraints['Embed-'+str(cnt_i)]

        #Create CNT again
        #THIS VERSION HAS THE OPTION profileNormal=ON
        CNT_Part(modelName, cnt_i, rad, acc_pts, acc_pts+N_p-1, cnt_coords)

        #Assign the CNT section to cnt_i
        mdb.models[modelName].parts[str_part].SectionAssignment(
            offset=0.0, 
            offsetField='', 
            offsetType=MIDDLE_SURFACE, 
            region=Region(cells=mdb.models[modelName].parts[str_part].cells),
            sectionName=cntMaterial, 
            thicknessAssignment=FROM_SECTION)

        #Generate CNT instance
        mdb.models[modelName].rootAssembly.Instance(
            dependent=ON, 
            name=str_inst,
            part=mdb.models[modelName].parts[str_part])
        
        plog('New part '+str_part+' and instance done')
        
        #Update previous cNT number
        cnt_prev = cnt_i

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

    #Hide the hollow box that was used to cut the GSs
    mdb.models[modelName].rootAssembly.features['CUTTER-1'].suppress()
    
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
    mdb.models[model].materials[materialName].Elastic(table=((elasModulus*1e9, poissonRatio), ))
    mdb.models[model].materials[materialName].Expansion(table=((expanCoefficient*1e-5, ), ))
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

#This function creates the CNT material and assigns it to a section
def Assign_Sections_CNTs(modelName, N_CNTs, cntSection):

    #Iterate over the CNTs
    for cnt_i in range(1, N_CNTs+1):

        #Get the string for the CNT part
        cnt_str = string_part('CNT', cnt_i)

        #Assign the CNT section to cnt_i
        mdb.models[modelName].parts[cnt_str].SectionAssignment(offset=0.0, offsetField='', offsetType=MIDDLE_SURFACE, 
	        region=Region(cells=mdb.models[modelName].parts[cnt_str].cells),
	        sectionName=cntSection, thicknessAssignment=FROM_SECTION)

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

#This function creates the assembly by creating the instances of each part
def Generate_CNT_Assembly(model, N_CNTs):

    #Create instances for each CNT
    for i in range(1, N_CNTs+1):

        #Generate name for CNT part
        str_cnt = string_part('CNT', i)

        #Generate name for CNT instance
        str_cnt_inst = string_instance('CNT', i)

        #Generate CNT instance
        mdb.models[model].rootAssembly.Instance(
            dependent=ON, name=str_cnt_inst,
            part=mdb.models[model].parts[str_cnt])	
 
#This function creates the sets that will be used when creating the embedded element constraints
def Sets_For_Embedded_Elements_CNTs(model, restartFlag, N_CNT_start, N_CNT_end, str_matrix, str_host):
    
    if restartFlag == 0:
        
        #Set for the matrix
        mdb.models[model].rootAssembly.Set(
    	    cells=mdb.models[model].rootAssembly.instances[str_matrix + '-1'].cells,
    	    name=str_host)

    #Sets for the CNTs
    for cnt_i in range(N_CNT_start, N_CNT_end+1):

        #Get the string for the CNT part
        cnt_str = string_part('CNT', cnt_i)

        #Generate name for CNT instance
        str_cnt_inst = string_instance('CNT', cnt_i)

        #Generate name of set for cnt_i
        set_str = ee_string('CNT', cnt_i)

        #Create the set for cnt_i
        mdb.models[model].rootAssembly.Set(
            cells=mdb.models[model].rootAssembly.instances[str_cnt_inst].cells,
            name=set_str)

#This functions creates the constraints for the embedded elements
def Embedded_Elements_Constraints_CNTs(model, N_CNT_start, N_CNT_end, str_matrix, str_host):

    #Iterate over the CNTs
    for cnt_i in range(N_CNT_start, N_CNT_end+1):

        #Get the string for the CNT set
        set_str = ee_string('CNT', cnt_i)

        #For cnt_i create an embedded element constraint
        #mdb.models['Model-1'].EmbeddedRegion(
        #	absoluteTolerance=0.0, fractionalTolerance=0.05, toleranceMethod=BOTH, weightFactorTolerance=1e-06,
        #	embeddedRegion=
        #		Region(cells=mdb.models['Model-1'].rootAssembly.instances[cnt_str + '-1'].cells.getSequenceFromMask(mask=('[#1 ]', ), )), 
        #	hostRegion=
        #		Region(cells=mdb.models['Model-1'].rootAssembly.instances[str_matrix + '-1'].cells.getSequenceFromMask(mask=('[#1 ]', ), )),
        #	name='EE-Constraint-%d' %(cnt_i))
        mdb.models[model].EmbeddedRegion(
        absoluteTolerance=0.0, 
        fractionalTolerance=0.05, 
        toleranceMethod=BOTH, 
        weightFactorTolerance=1e-06,
        embeddedRegion=mdb.models[model].rootAssembly.sets[set_str],
        hostRegion=mdb.models[model].rootAssembly.sets[str_host],
        name='EE-CNT-%d' %(cnt_i))
 
#This function creates the sets that will be used when creating the embedded element constraints
def Sets_For_Embedded_Elements_Array_CNTs(model, cnts_new):

    #Sets for the CNTs
    for cnt_i in cnts_new:

        #Get the string for the CNT part
        cnt_str = string_part('CNT', cnt_i)

        #Generate name for CNT instance
        str_cnt_inst = string_instance('CNT', cnt_i)

        #Generate name of set for cnt_i
        set_str = ee_string('CNT', cnt_i)

        #Create the set for cnt_i
        mdb.models[model].rootAssembly.Set(
            cells=mdb.models[model].rootAssembly.instances[str_cnt_inst].cells,
            name=set_str)

#This functions creates the constraints for the embedded elements
def Embedded_Elements_Constraints_Array_CNTs(model, cnts_new, str_matrix, str_host):

    #Iterate over the CNTs
    for cnt_i in cnts_new:

        #Get the string for the CNT set
        set_str = ee_string('CNT', cnt_i)

        #For cnt_i create an embedded element constraint
        mdb.models[model].EmbeddedRegion(
            absoluteTolerance=0.0, 
            fractionalTolerance=0.05, 
            toleranceMethod=BOTH, 
            weightFactorTolerance=1e-06,
            embeddedRegion=mdb.models[model].rootAssembly.sets[set_str],
            hostRegion=mdb.models[model].rootAssembly.sets[str_host],
            name='EE-CNT-%d' %(cnt_i))

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
def Create_Embedded_Elements_GSs(model, matrixName, gs_i, gs_inst_str):

    #Add embedded region to GS i
    mdb.models[model].EmbeddedRegion(
        absoluteTolerance=0.0, 
        embeddedRegion=Region(cells=mdb.models[model].rootAssembly.instances[gs_inst_str].cells),
        fractionalTolerance=0.05, 
        hostRegion=Region(cells=mdb.models[model].rootAssembly.instances[matrixName + '-1'].cells),
        name='EE-GS-%d' %(gs_i),
        toleranceMethod=BOTH, 
        weightFactorTolerance=1e-06)

#This function generates all meshes
#Embedded element constraints for GSs are added here because GS may be trimmed 
def Generate_Meshes(modelName, matrixName, selectedElementCode, eeMeshSize, numGSs, indexOutside, N_CNTs, cnt_struct, cnt_coords):
    
    #Generate the meshes for the GSs
    Mesh_all_GSs(modelName, matrixName, selectedElementCode, eeMeshSize, numGSs, indexOutside)

    #Generate meshes for the CNTs
    Generate_CNT_Meshes(modelName, N_CNTs, cnt_struct, cnt_coords)

    #Generate the mesh for the matrix
    Generate_Matrix_Mesh(modelName, matrixName, selectedElementCode)
        
    #This seems to be required by Abaqus
    mdb.models[modelName].rootAssembly.regenerate()

#Mesh all GSs
def Mesh_all_GSs(modelName, matrixName, selectedElementCode, eeMeshSize, numGSs, indexOutside):
    
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
        Create_Embedded_Elements_GSs(modelName, matrixName, gs_i, gs_inst_str)
        
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

#Mesh single GS
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

#This function generates the mesh for the CNTs
def Generate_CNT_Meshes(modelName, N_CNTs, cnt_struct, cnt_coords):

    #Number of accumulated points
    acc_pts = 0

    #Go through every CNT to mesh each of them
    for cnt_i in range(1, N_CNTs+1):

        #Number of points in CNTi and its radius
        N_p = int(cnt_struct[cnt_i][0])
        cnt_rad = cnt_struct[cnt_i][1]

        #Get the string for the CNT part
        cnt_str = string_part('CNT', cnt_i)
        
        #Size of element
        cnt_el = 2.0*cnt_struct[cnt_i][1]

        #Mesh cnt_i, use its radius as the element size
        #deviationFactor and minSizeFactor have the default values from Abaqus
        mdb.models[modelName].parts[cnt_str].seedPart(
            deviationFactor=0.1, minSizeFactor=0.1, size=cnt_el)
        mdb.models[modelName].parts[cnt_str].generateMesh()

        #Check if the CNT needs to be remeshed
        Remesh_CNT_When_Needed(modelName, cnt_str, cnt_i, cnt_rad, acc_pts, N_p, cnt_el, cnt_coords)

        #Increase the number of accumulated points
        acc_pts += N_p

#This function checks for the cases where a CNT needs to be remeshed
def Remesh_CNT_When_Needed(modelName, cnt_str, cnt_i, cnt_rad, acc_pts, N_p, cnt_el, cnt_coords):

    #Get the number of nodes generated
    num_nodes = mdb.models[modelName].parts[cnt_str].getMeshStats((mdb.models[modelName].parts[cnt_str].cells,)).numNodes
    #Get the number of hexahedral elements generated
    num_els = mdb.models[modelName].parts[cnt_str].getMeshStats((mdb.models[modelName].parts[cnt_str].cells,)).numHexElems
    #Calculate the minimum theoretical number of elements
    #The theoretical number of elements is 4(N_p-1), since there are 4 elements for each CNT segment (a CNT has N_p-1 segments)
    theor_els = 4*(N_p-1)
    
    #Check the number of nodes in the mesh
    if num_nodes == 0:

        #The CNT was not meshed
        #Cut the CNT cell where needed
        Partition_CNT_Cell(modelName, cnt_rad, acc_pts, acc_pts+N_p-1, cnt_coords, cnt_str)

        #Try to mesh again
        mdb.models[modelName].parts[cnt_str].seedPart(deviationFactor=0.1, minSizeFactor=0.1, size=cnt_el)
        mdb.models[modelName].parts[cnt_str].generateMesh()
        
    #Check if number of elements is within limits, if not within limits probably ppart of the CNT is a beam section
    #This beacause some CNT segments are longer close to the boundaries because a short segment was merged with a regular segment
    #In such cases Abaqus creates 8 elements in that segment instead of only 4
    #Thus, if a CNT is cut by the boundary by one side it may have 4 elements more and if it is cut twice it could have 
    #up to 8 elements more
    elif num_els < theor_els:
        
        #Delete part and re-do it with the flag for profile normal constant
        plog('{} nodes={} num_els={} 4*(N_p-1)={}\n'.format(cnt_str,num_nodes,num_els,theor_els))
        Delete_And_CreateAgain(modelName, cnt_str, cnt_i, cnt_rad, acc_pts, acc_pts+N_p-1, cnt_coords)

        #Try to mesh again
        mdb.models[modelName].parts[cnt_str].seedPart(deviationFactor=0.1, minSizeFactor=0.1, size=cnt_el)
        mdb.models[modelName].parts[cnt_str].generateMesh()

#This function cuts a CNT cell when it could not be meshed
def Partition_CNT_Cell(modelName, cnt_rad, cnt_start, cnt_end, cnt_coords, str_part):

    #Half of the cylinder height
    hc = cnt_rad*0.1
    
    #Array to store additional cuts
    additionalCuts = []
    
    #Array to store cuts
    cutsDone = []
    
    #Array to store all cuts that could not be done at first
    invalidCuts = []

    #Iterate over the CNT points, starting on the secont point and finishing on the previous to last point
    for i in range(cnt_start+1, cnt_end): 

        #print('i=',i-cnt_start,' cells=', len(mdb.models[modelName].parts[str_part].cells))

        #Get the points of the CNT segments that share point i
        #End point of first CNT segment
        P1 = cnt_coords[i-1]
        #Point shared by both CNT segments
        P2 = cnt_coords[i]
        #End point of second CNT segment
        P3 = cnt_coords[i+1]

        #Get the unit vector of first CNT segment
        v1 = get_unit_vector(P2, P1)

        #Get the unit vector of second CNT segment
        v2 = get_unit_vector(P2, P3)

        #Calculate the dot product of v1 and v2 to obtain the cosine of the angle between them
        #Need the dot product with negative sign
        cosV = - dot(v1,v2)

        #Check if need to cut the cell at point P2
        if cosV < cos45:
            #The CNT cell needs to be cut at P2
        
            #Add index to array of cuts done
            cutsDone.append(i)
            
            #Check if the current point was added to the array of additional cuts
            if len(additionalCuts) > 0 and i == additionalCuts[-1]:
                #Remove the last element from additionalCuts array
                additionalCuts.pop()
                
            #Get the edges of the octagon that are used to cut the CNT
            octagon_edges = Get_octagon_to_cut_CNT(modelName, str_part, cnt_rad, P2, v1, v2, hc)

            #Datum points for testing
            #if i == cnt_start+1:
            #    for j in range(len(octagon)):
            #        mdb.models[modelName].parts[str_part].DatumPointByCoordinate(coords=octagon[j].pointOn[0])
            #Use the octagon edges in the tuple to partition the cell
            # mdb.models[modelName].parts[str_part].PartitionCellByPatchNEdges(
            #     cell=mdb.models[modelName].parts[str_part].cells.findAt(P3, ),
            #     edges=octagon_edges)
            
            try:
                #Use the octagon edges in the tuple to partition the cell
                mdb.models[modelName].parts[str_part].PartitionCellByPatchNEdges(
                    cell=mdb.models[modelName].parts[str_part].cells.findAt(P3, ),
                    edges=octagon_edges)
            except:
                #Datum point could be helpful for debugging
                #mdb.models[modelName].parts[str_part].DatumPointByCoordinate(P2)
                
                #Save the point for invalid cut
                invalidCuts.append(i)
                
                #Additional couts are needed, so append the previous and next point number
                #to the array of additional cuts
                #Check the previous point is not the first CNT point, which happens when:
                #   the point is not the first point in the CNT
                #   AND
                #   the poit is not a previus cut (i.e., cutsDone is empty or point is not the last element in array cutsDone)
                if i-1 != cnt_start and (len(cutsDone)==0 or i-1 != cutsDone[-1]):
                    additionalCuts.append(i-1)
                #Check the next point is not the last point in the CNT
                if i+1 < cnt_end-1:
                    additionalCuts.append(i+1)
                
                plog('Could not cut part {} on point {}, cnt_start={} cnt_end={}\n'.format(str_part, i-cnt_start, cnt_start, cnt_end))
                
    #Check if there are additional cuts needed
    if len(additionalCuts) > 0:
    
        #Addtional cuts needed, iterate over those points
        for j in additionalCuts:
    
            #Get the points of the CNT segments that share point i
            #End point of first CNT segment
            P1 = cnt_coords[j-1]
            #Point shared by both CNT segments
            P2 = cnt_coords[j]
            #End point of second CNT segment
            P3 = cnt_coords[j+1]
    
            #Get the unit vector of first CNT segment
            v1 = get_unit_vector(P2, P1)
    
            #Get the unit vector of second CNT segment
            v2 = get_unit_vector(P2, P3)
    
            #Get the edges of the octagon that are used to cut the CNT
            octagon_edges = Get_octagon_to_cut_CNT(modelName, str_part, cnt_rad, P2, v1, v2, hc)
            try:
                #Use the octagon edges in the tuple to partition the cell
                mdb.models[modelName].parts[str_part].PartitionCellByPatchNEdges(
                    cell=mdb.models[modelName].parts[str_part].cells.findAt(P2, ),
                    edges=octagon_edges)
    
                plog('Additional cut on part {} on point {}, cnt_start={} cnt_end={}\n'.format(str_part, j-cnt_start, cnt_start, cnt_end))
            except:
    
                plog('Could not make additional cut on part {} on point {}, cnt_start={} cnt_end={}\n'.format(str_part, j-cnt_start, cnt_start, cnt_end))
    
    #Check if the invalid cuts can be done now
    if len(invalidCuts) > 0:  
    
        #Addtional cuts needed, iterate over those points
        for j in invalidCuts:
            #print('Invalid cut second attempt ',j)
    
            #Get the points of the CNT segments that share point i
            #End point of first CNT segment
            P1 = cnt_coords[j-1]
            #Point shared by both CNT segments
            P2 = cnt_coords[j]
            #End point of second CNT segment
            P3 = cnt_coords[j+1]
    
            #Get the unit vector of first CNT segment
            v1 = get_unit_vector(P2, P1)
    
            #Get the unit vector of second CNT segment
            v2 = get_unit_vector(P2, P3)
    
            #Get the edges of the octagon that are used to cut the CNT
            octagon_edges = Get_octagon_to_cut_CNT(modelName, str_part, cnt_rad, P2, v1, v2, hc)
            try:
                #Use the octagon edges in the tuple to partition the cell
                mdb.models[modelName].parts[str_part].PartitionCellByPatchNEdges(
                    cell=mdb.models[modelName].parts[str_part].cells.findAt(P2, ),
                    edges=octagon_edges)
    
                plog('Cut after second attempt on part {} on point {}, cnt_start={} cnt_end={}\n'.format(str_part, j-cnt_start, cnt_start, cnt_end))
            except:
    
                plog('Could not make cut after second attempt on part {} on point {}, cnt_start={} cnt_end={}\n'.format(str_part, j-cnt_start, cnt_start, cnt_end))     

#This function get the octagon needed to cut a CNT
def Get_octagon_to_cut_CNT(modelName, str_part, cnt_rad, P2, v1, v2, hc):
    
    #Get a unit vector that goes along the plane that bisects the angle 
    #between v1 and v2
    vm = make_unit((v1[0]+v2[0], v1[1]+v2[1], v1[2]+v2[2]))
    #print('vm')

    #Vector normal to v1, v2, and vm (all three lie on the same plane)
    N3 = cross(v2, v1)

    #Get the normal for the plane at the midpoint
    #Since both v1 and vm are unit vectors, N is also a unit vector
    N = cross(N3, vm)

    #Calculate a displacement to set the points of the cylinder
    disp = (N[0]*hc, N[1]*hc, N[2]*hc)

    #Calculate first point of Cylinder height
    C1 = (P2[0]+disp[0], P2[1]+disp[1], P2[2]+disp[2])

    #Calculate second point of Cylinder height
    C2 = (P2[0]-disp[0], P2[1]-disp[1], P2[2]-disp[2])

    #Calculate vector along the plane
    #Since both v1 and vm are unit vectors, S is also a unit vector
    S = cross(N3, v1)

    #Calculate the radius of the cylinder that will be used to select the points needed to cut the cell
    #Also increase the radius 0.5% so that vertices are not missed die to numerical erros
    rad_cyl = cnt_rad*1.005/(-dot(vm, S))
    #print('rad_cyl=',rad_cyl)

    #Select the edges enclosed by the cylinder with centerpoints C1 and C2 and radius rad_cyl
    octagon = mdb.models[modelName].parts[str_part].edges.getByBoundingCylinder(center1=C1, center2=C2, radius=rad_cyl)
    #print('i=',i-cnt_start-1," octagon.len=",len(octagon))
    #print(octagon)
    #print(octagon[0])
    #print(octagon[1])
    #For some reason, the selected edges are not a tuple
    #Thus, create a tuple with the edges of the octagon
    octagon_edges = (octagon[0], octagon[1], octagon[2], octagon[3], octagon[4], octagon[5], octagon[6], octagon[7])
    
    return octagon_edges

#This function deletes a CNt that could not be meshed beacause part of it was generated as a shell element
def Delete_And_CreateAgain(modelName, str_part, cnt_i, cnt_rad, cnt_start, cnt_end, cnt_coords):
    
    #Delete old embedded element constraint
    del mdb.models[modelName].constraints['EE-CNT-%d' %(cnt_i)]

    #Generate name of set for embedded element constraint for cnt_i
    set_str = ee_string('CNT', cnt_i)

    #Delete old set for cnt_i
    del mdb.models[modelName].rootAssembly.sets[set_str]

    #Generate name for CNT instance
    str_inst = string_instance('CNT', cnt_i)
    
    #Delete old instance
    del mdb.models[modelName].rootAssembly.features[str_inst]
    
    #Delete initial sweep
    del mdb.models[modelName].parts[str_part].features['Solid sweep-1']
    
    plog('Sweep feature of part '+str_part+' and instance deleted\n')
    plog('cnt_end= {} cnt_start={}\n'.format(cnt_end, cnt_start))

    #Create sweep for CNT again

    #Profile name 
    profile_str = '__profile__'
    
    #Generate the octagon that will be swept through the CNT segments
    Generate_Sweep_Profile(modelName, cnt_coords[cnt_end-1], cnt_coords[cnt_end], cnt_rad, str_part, profile_str)
    
    #Perform Sweep  
    Sweep_Command_profileNormal_ON(modelName, str_part, profile_str, cnt_start, cnt_end)

    #Delete sketch
    del mdb.models[modelName].sketches[profile_str]

    #Assign the CNT section to cnt_i
    mdb.models[modelName].parts[str_part].SectionAssignment(
        offset=0.0, 
        offsetField='', 
        offsetType=MIDDLE_SURFACE, 
        region=Region(cells=mdb.models[modelName].parts[str_part].cells),
        sectionName=cntMaterial, 
        thicknessAssignment=FROM_SECTION)

    #Generate CNT instance
    mdb.models[modelName].rootAssembly.Instance(
        dependent=ON, 
        name=str_inst,
        part=mdb.models[modelName].parts[str_part])

    #Create new set for cnt_i
    mdb.models[modelName].rootAssembly.Set(
        cells=mdb.models[modelName].rootAssembly.instances[str_inst].cells,
        name=set_str)
    
    #Create new embedded element constraint
    mdb.models[modelName].EmbeddedRegion(
        absoluteTolerance=0.0, 
        fractionalTolerance=0.05, 
        toleranceMethod=BOTH, 
        weightFactorTolerance=1e-06,
        embeddedRegion=mdb.models[modelName].rootAssembly.sets[set_str],
        hostRegion=mdb.models[modelName].rootAssembly.sets[strHost],
        name='EE-CNT-%d' %(cnt_i))
    
    plog('New part '+str_part+' and instance done\n')

#This function generates the mesh for the CNTs
def Generate_CNT_Meshes_Range(modelName, N_CNT_start, N_CNT_end, cnt_struct, cnt_coords):
    
    #Number of accumulated points
    acc_pts = Accumulate_points(N_CNT_start, cnt_struct)

    #Go through every CNT to mesh each of them
    for cnt_i in range(N_CNT_start, N_CNT_end+1):

        #Number of points in CNTi and its radius
        N_p = int(cnt_struct[cnt_i][0])
        cnt_rad = cnt_struct[cnt_i][1]

        #Get the string for the CNT part
        cnt_str = string_part('CNT', cnt_i)
        
        #Size of element
        cnt_el = 2.0*cnt_struct[cnt_i][1]

        #Mesh cnt_i, use its radius as the element size
        #deviationFactor and minSizeFactor have the default values from Abaqus
        mdb.models[modelName].parts[cnt_str].seedPart(
            deviationFactor=0.1, minSizeFactor=0.1, size=cnt_struct[cnt_i][1])
        mdb.models[modelName].parts[cnt_str].generateMesh()

        #Check if the CNT needs to be remeshed
        Remesh_CNT_When_Needed(modelName, cnt_str, cnt_i, cnt_rad, acc_pts, N_p, cnt_el, cnt_coords)

        #Increase the number of accumulated points
        acc_pts += N_p
        
    #This seems to be required by Abaqus
    mdb.models[modelName].rootAssembly.regenerate()

#This function generates the mesh for the CNTs
def Generate_CNT_Meshes_Array(modelName, cnts_new, arr_acc_pts, cnt_struct, cnt_coords):
    
    #Number of CNTs to mesh
    n_cnts = len(cnts_new)

    #Go through every CNT to mesh each of them
    for i in range(n_cnts):
        
        #Get CNT number
        cnt_i = cnts_new[i]

        #Number of points in CNTi and its radius
        N_p = int(cnt_struct[cnt_i][0])
        cnt_rad = cnt_struct[cnt_i][1]

        #Get the string for the CNT part
        cnt_str = string_part('CNT', cnt_i)
        
        #Size of element
        cnt_el = 2.0*cnt_struct[cnt_i][1]

        #Mesh cnt_i, use its radius as the element size
        #deviationFactor and minSizeFactor have the default values from Abaqus
        mdb.models[modelName].parts[cnt_str].seedPart(
            deviationFactor=0.1, minSizeFactor=0.1, size=cnt_struct[cnt_i][1])
        mdb.models[modelName].parts[cnt_str].generateMesh()

        #Check if the CNT needs to be remeshed
        Remesh_CNT_When_Needed(modelName, cnt_str, cnt_i, cnt_rad, arr_acc_pts[i], N_p, cnt_el, cnt_coords)
        
    #This seems to be required by Abaqus
    mdb.models[modelName].rootAssembly.regenerate()

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

#This function creates all sets for the nodes that correspond to the centerline of that CNT
def Create_All_Sets_For_CNT_Points(modelName, N_CNTs, cnt_struct, cnt_coords):

    #Initializae the number of accumulated points to zero
    acc_pts = 0

    #Iterate over all CNTs
    for cnt_i in range(1, N_CNTs+1):

        #Number of points in CNTi and its radius
        N_p = int(cnt_struct[cnt_i][0])
        cnt_rad = cnt_struct[cnt_i][1]

        #Get the string for part corresponding to cnt_i
        cnt_str = string_instance('CNT', cnt_i)

        #Create a set for the nodes at the centerline of cnt_i
        Create_Set_For_CNT_Points(modelName, cnt_i, cnt_rad, acc_pts, acc_pts+N_p-1, cnt_coords, cnt_str)

        #Increase the number of accumulated points
        acc_pts += N_p

#This function creates a set for the nodes that correspond to the centerline of that CNT
def Create_Set_For_CNT_Points(modelName, cnt_i, cnt_rad, cnt_start, cnt_end, cnt_coords, cnt_str):

    #Calculate radius for serching nodes
    new_rad = cnt_rad/10

    #Get the name of the node set
    node_set_str = cnt_string_node_set(cnt_i)
    
    #Array to store all nodes that correspond to CNT points
    allNodes = []
    
    #Array to store the labels of nodes when more than one is founnd using getByBoundingSphere
    labelsToRemove = []

    #Iterate over all nodes in the centerline of the CNT
    for i in range(cnt_start, cnt_end+1):
        
        #Find CNT node (or sometimes nodes)
        newNodes = mdb.models[modelName].rootAssembly.instances[cnt_str].nodes.getByBoundingSphere(
            center=cnt_coords[i],
            radius=new_rad)
        
        #Append node to array
        allNodes.append(newNodes)
        
        #Check if more than one node was found (sometimes this happens)
        if len(newNodes) > 1:
            
            #Get the number of extra nodes
            extraNodes = len(newNodes) - 1
            
            #Append the labels of the extra nodes to the array labelsToRemove
            for j in range(extraNodes):
                #Ignore the first node, so the +1 ensures iterating from 1
                labelsToRemove.append(newNodes[j+1].label)
    
    #Create a set with the name of the node set and containing all nodes in the centerline of the CNT
    mdb.models[modelName].rootAssembly.Set(
        nodes=allNodes,
        name=node_set_str)

    #Print the length of the set
    #plog('%s nodes=%d points=%d'%(node_set_str, len(mdb.models[modelName].rootAssembly.sets[node_set_str].nodes), cnt_end+1-cnt_start))
    
    #Check if there are nodes to remove
    if len(labelsToRemove) >= 1:
        
        #Create a temporary set with the node labels that need to be removed
        setToRemove = mdb.models[modelName].rootAssembly.SetFromNodeLabels(
            name='TMP', 
            nodeLabels=( (cnt_str, labelsToRemove), ))
        
        #Remove them from the original set
        mdb.models[modelName].rootAssembly.SetByBoolean(
            name=node_set_str, 
            sets=(cntSet, setToRemove, ), 
            operation=DIFFERENCE)
        
        #Delete temporary set
        del mdb.models[modelName].rootAssembly.sets['TMP']
        
        #Print the length of the set
        #plog('%s nodes=%d points=%d'%(node_set_str, len(mdb.models[modelName].rootAssembly.sets[node_set_str].nodes), cnt_end+1-cnt_start))

#This function creates all sets for the nodes that correspond to the centerline of that CNT
def Create_All_Sets_For_CNT_Points_Range(modelName, N_CNT_start, N_CNT_end, cnt_struct, cnt_coords):
    
    #Number of accumulated points
    acc_pts = Accumulate_points(N_CNT_start, cnt_struct)

    #Iterate over all CNTs
    for cnt_i in range(N_CNT_start, N_CNT_end+1):

        #Number of points in CNTi and its radius
        N_p = int(cnt_struct[cnt_i][0])
        cnt_rad = cnt_struct[cnt_i][1]

        #Get the string for part corresponding to cnt_i
        cnt_str = string_instance('CNT', cnt_i)

        #Create a set for the nodes at the centerline of cnt_i
        Create_Set_For_CNT_Points(modelName, cnt_i, cnt_rad, acc_pts, acc_pts+N_p-1, cnt_coords, cnt_str)

        #Increase the number of accumulated points
        acc_pts += N_p 

#This function creates all sets for the nodes that correspond to the centerline of that CNT
def Create_All_Sets_For_CNT_Points_Array(modelName, cnts_new, arr_acc_pts, cnt_struct, cnt_coords):
    
    #Get the number of CNTs
    n_cnts = len(cnts_new)

    #Iterate over all CNTs
    for i in range(n_cnts):
        
        #Get the CNT number
        cnt_i = cnts_new[i]

        #Number of points in CNTi and its radius
        N_p = int(cnt_struct[cnt_i][0])
        cnt_rad = cnt_struct[cnt_i][1]

        #Get the string for part corresponding to cnt_i
        cnt_str = string_instance('CNT', cnt_i)
        
        #Get the number of accumulated points
        acc_pts = arr_acc_pts[i]

        #Create a set for the nodes at the centerline of cnt_i
        Create_Set_For_CNT_Points(modelName, cnt_i, cnt_rad, acc_pts, acc_pts+N_p-1, cnt_coords, cnt_str)

#This function creates an element set that contains the elements in the extended region of the RVE
#that need to be hidden in the visualization
def Sets_For_Elements_To_Hide(modelName, matrixName, P0, Lxyz, matrixMeshSize, halfMatrixMeshSize):
    
    #Initialize empty array
    elsToHide = []
    
    #Name for the set that will contain all elements in the extended layer
    hideSetName = 'HIDE-SET'
    
    #String for matrix instance
    matrixInstance = matrixName + '-1'
    
    #Add elements from top
    elsToHide.append(mdb.models[modelName].rootAssembly.instances[matrixInstance].elements.getByBoundingBox(
        P0[0]-1.5*matrixMeshSize, P0[1]-1.5*matrixMeshSize, P0[2]+Lxyz[2]-halfMatrixMeshSize,
        P0[0]+Lxyz[0]+1.5*matrixMeshSize, P0[1]+Lxyz[1]+1.5*matrixMeshSize, P0[2]+Lxyz[2]+1.5*matrixMeshSize))
    #print('len(elsToHide)=',len(elsToHide))
    
    #Add elements from bottom
    elsToHide.append(mdb.models[modelName].rootAssembly.instances[matrixInstance].elements.getByBoundingBox(
        P0[0]-1.5*matrixMeshSize, P0[1]-1.5*matrixMeshSize, P0[2]-1.5*matrixMeshSize,
        P0[0]+Lxyz[0]+1.5*matrixMeshSize, P0[1]+Lxyz[1]+1.5*matrixMeshSize, P0[2]+halfMatrixMeshSize))
    #print('len(elsToHide)=',len(elsToHide))
    
    #Add elements from sides
    #Height along y-axis
    elsToHide.append(mdb.models[modelName].rootAssembly.instances[matrixInstance].elements.getByBoundingBox(
        P0[0]-1.5*matrixMeshSize, P0[1]+Lxyz[1]-halfMatrixMeshSize, P0[2]-1.5*matrixMeshSize,
        P0[0]+Lxyz[0]+1.5*matrixMeshSize, P0[1]+Lxyz[1]+1.5*matrixMeshSize, P0[2]+Lxyz[2]+1.5*matrixMeshSize))
    #print('len(elsToHide)=',len(elsToHide))
    elsToHide.append(mdb.models[modelName].rootAssembly.instances[matrixInstance].elements.getByBoundingBox(
        P0[0]-1.5*matrixMeshSize, P0[1]-1.5*matrixMeshSize, P0[2]-1.5*matrixMeshSize,
        P0[0]+Lxyz[0]+1.5*matrixMeshSize, P0[1]+halfMatrixMeshSize, P0[2]+Lxyz[2]+1.5*matrixMeshSize))
    #print('len(elsToHide)=',len(elsToHide))
    #Height along x-axis
    elsToHide.append(mdb.models[modelName].rootAssembly.instances[matrixInstance].elements.getByBoundingBox(
        P0[0]+Lxyz[0]-halfMatrixMeshSize, P0[1]-1.5*matrixMeshSize, P0[2]-1.5*matrixMeshSize,
        P0[0]+Lxyz[0]+1.5*matrixMeshSize, P0[1]+Lxyz[1]+1.5*matrixMeshSize, P0[2]+Lxyz[2]+1.5*matrixMeshSize))
    #print('len(elsToHide)=',len(elsToHide))
    elsToHide.append(mdb.models[modelName].rootAssembly.instances[matrixInstance].elements.getByBoundingBox(
        P0[0]-1.5*matrixMeshSize, P0[1]-1.5*matrixMeshSize, P0[2]-1.5*matrixMeshSize,
        P0[0]+halfMatrixMeshSize, P0[1]+Lxyz[1]+1.5*matrixMeshSize, P0[2]+Lxyz[2]+1.5*matrixMeshSize))
    #print('len(elsToHide)=',len(elsToHide))
    
    mdb.models[modelName].rootAssembly.Set(elements=elsToHide, name=hideSetName)

#Sets needed to create the PBC equations
def Create_Set_for_PBC(model, matrixName, P0, Lxyz, corner):

    modelRoot = mdb.models[model].rootAssembly
    
    #Calculate half the lengths of the RVE along each direction
    half_x = Lxyz[0]*0.5
    half_y = Lxyz[1]*0.5
    half_z = Lxyz[2]*0.5
    
    #Set containing the 6 faces of the RVE (temperature will be applied to this set)
    #print('External_Faces')
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
    #print('Individual faces')
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


##############################################################################################
##############################################################################################
##############################################################################################
#####################################---MAIN PROGRAM---#######################################


######################################---OPENING CSV---#######################################

start0 = time.time()

#Read CNT point coordinates
with open(coord_file) as f:
    #x = csv.reader(f, quoting=csv.QUOTE_NONNUMERIC)
    #cnt_coords = list(x)
    cnt_coords = list(csv.reader(f, quoting=csv.QUOTE_NONNUMERIC))
#print(cnt_coords[0])
#print(cnt_coords[0][0])

#Read the number of points per CNT and the radii
with open(struct_file) as f:
    #x = csv.reader(f, quoting=csv.QUOTE_NONNUMERIC)
    #cnt_struct = list(x)
    cnt_struct = list(csv.reader(f, quoting=csv.QUOTE_NONNUMERIC))
#print(cnt_struct[0])
#print(int(cnt_struct[0][0]))

#Get the number of CNTs
N_CNTs = int(cnt_struct[0][0])

#Array of CNTs to delete and re-do
cnts_new = [402]

#Name of the job to be used based on its parameters
jobName = 'VOL0160_CNT01GNP02'

#Name of the file to save print messages
print_file = jobName + '.txt'
pfile = open(print_file, "a")

plog('New input file: '+jobName+'.inp\n')
    
start = time.time()

#Import model from inp
inpName = jobName+'.inp'
plog('Reading input file: ' + inpName + ' ...\n')
mdb.ModelFromInputFile(
    inputFileName=inpName, 
    name=modelName)

plog('Importing model from inp time: {} secs.\n'.format(time.time()-start))

#Select an element type based on the choice of geometric order
selectedElementCode = Select_Elemet_Type(selectedGeomOrder)

####################################---MODEL CREATION---######################################

#Creating a section for each material
Create_Section(modelName, matrixMaterial)
Create_Section(modelName, gsMaterial)
Create_Section(modelName, cntMaterial)

#-----------------------------------START: CREATE ASSEMBLY-----------------------------------#

start = time.time()

#This line seems to be required by Abaqus to set a Cartesian system
mdb.models[modelName].rootAssembly.DatumCsysByDefault(CARTESIAN)

#Array to store accumulated points
arr_acc_pts = []

#Delete and re-do the partts in array cnt_new
Delete_Imported_Create_New(modelName, cnts_new, cnt_struct, cnt_coords, cntMaterial, arr_acc_pts)

end = time.time()
plog("Time for part and instance generation: {} secs.\n".format(end-start))

#------------------------------------END: CREATE ASSEMBLY------------------------------------#

#-----------------------------START: PERIODIC BOUNDARY CONDITIONS----------------------------#

start = time.time()

#Only create the sets for new CNTs

#Create sets that will be used when creating the embedded element constraints for CNTs
Sets_For_Embedded_Elements_Array_CNTs(modelName, cnts_new)

#Create embedded elements constraints for CNTs
Embedded_Elements_Constraints_Array_CNTs(modelName, cnts_new, matrixName, strHost)

end = time.time()
plog("Time for embedded sets and constraints: {} secs.\n".format(end-start))

#------------------------------END: PERIODIC BOUNDARY CONDITIONS-----------------------------#

#---------------------------------------START: MESHING---------------------------------------#

start = time.time()
    
#Mesh only the new CNTs in the range
Generate_CNT_Meshes_Array(modelName, cnts_new, arr_acc_pts, cnt_struct, cnt_coords)

#Create sets for central CNT nodes only the new CNTs in the range
#NOTE: Sets are generated on root assembly
Create_All_Sets_For_CNT_Points_Array(modelName, cnts_new, arr_acc_pts, cnt_struct, cnt_coords)

end = time.time()
plog("Time for meshing: {} secs.\n".format(end-start))

#print('The model has been completed in %s seconds.' % round(time.time() - start0, 1))
#---------------------------------------END: MESHING----------------------------------------#

###################################---JOB CREATION---######################################

    
#Create and submit job using Abaqus default values
mdb.Job(atTime=None, contactPrint=OFF, description='', echoPrint=OFF, 
    explicitPrecision=SINGLE, getMemoryFromAnalysis=True, historyPrint=OFF, 
    memory=90, memoryUnits=PERCENTAGE, model=modelName, modelPrint=OFF, 
    multiprocessingMode=DEFAULT, name=jobName, nodalOutputPrecision=SINGLE, 
    numCpus=int(multiprocessing.cpu_count()*0.8), numDomains=int(multiprocessing.cpu_count()*0.8), numGPUs=1, queue=None, resultsFormat=ODB, scratch=
    '', type=ANALYSIS, userSubroutine='', waitHours=0, waitMinutes=0)
        
#Save input file
start = time.time()
mdb.jobs[jobName].writeInput(consistencyChecking=OFF)
plog("Input file written: {} secs.\n".format(time.time()-start))
     
plog("Time for Abaqus model: {} secs.\n".format(time.time() -start0))
pfile.close()