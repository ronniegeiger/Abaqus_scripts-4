# -*- coding: mbcs -*-
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

#Pyton libraries
import math
import csv
import time
import sys

######################################---GLOBAL VARIABLES---########################################

#Files to read
coord_file = 'cnt_coordinates.csv'
struct_file = 'cnt_struct.csv'
sample_file = 'sample_geom.csv'

###################### Mesh flags and variables
#Number of elements per side of the RVE
elementsPerSide = 10
#EmbeddedMesh/HostMesh ratio
meshRatio = 0.5
    
#Define the string for the name of the matrix part
matrixName = 'MATRIX'

#String for the host set (i.e., the matrix)
str_host = 'HOST_SET'

#Cosine of 45 degrees (PI/4)
cos45 = 0.7071067812

#"Zero" for comparing floating point numbers
Zero = 1e-7

#Displacement to be applied (microns)
disp = 0.5

#Matrix properties
#Define the name and section of the matrix
matrixMaterial = 'MATRIX_MAT'
matrixSection = 'MATRIX_SEC'
#Define the mass density (kg/m3)
matrixDensity = 905
#Define the elastic modulus (MPa)
matrixModulus = 2000
#Define the Poisson ratio
matrixPoissonR = 0.42
#Define de coefficient of thermal expansion (e-5 1/C)
matrixExpCoeff = 10.5
#Define the thermal conductivity (W/m*K)
matrixThermConductivity = 0.19
#Define the specific heat (J/mol*K)
matrixSpecHeat = 75
#Define the electrical conductivity (S/m)
matrixElecConductivity = 200

#CNT properties
#Define the name and section of the filler
cntMaterial = 'CNT_MAT'
cntSection = 'CNT_SEC'
#Define the mass density (kg/m3)
cntDensity = 2200
#Define the elastic modulus (MPa)
cntModulus = 1e6
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
cnt_rad_max = 0.03

#Some names
modelName = 'Model-1'
stpName = 'Step-1'

#Increment for step
inicialIncrement=0.1

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
submitJob = 0

#Displacement flags
#These flags inidicate if displacement is applied in a direction
# False = no displacement
# True = displacement as indicated in the variable for displacement
dispXflag = True
dispYflag = False
dispZflag= False
dispX = 0.15
dispY = 0.0
dispZ = 0.0

######################################---GLOBAL VARIABLES---########################################
######################################----------END---------########################################

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

#This function generates a string for the CNT part
def cnt_string_part(cnt_i):
	#Return the string in the format CNT-cnt_i
	return 'CNT-%d' %(cnt_i)

#This function generates a string for the CNT instance
def cnt_string_instance(cnt_i):
	#Return the string in the format CNT-cnt_i
	return 'CNT-%d-1' %(cnt_i)

#This function generates a string for the node set of CNT i
def cnt_string_node_set(cnt_i):
	#Return the string in the format CNT-nodes-cnt_i
	return 'CNT-%d-nodes' %(cnt_i)

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

######################################---ABAQUS FUNCTIONS---########################################

def matrix_part(modelName, P0, Lxyz, Lxyz_ext, matrixName):
	
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
	
#This function generates all edges of a CNT
def generate_edges(modelName, cnt_start, cnt_end, cnt_coords, str_part):
	
	#Create edges, iterate over all points of the CNT
	for i in range(cnt_start, cnt_end):
	#for i in range(cnt_start, cnt_start+11):
		
		#First point of the CNT segment
		P1 = cnt_coords[i]
		#Second point of the CNT segment
		P2 = cnt_coords[i+1]
		
		#Generate the segment
		mdb.models[modelName].parts[str_part].WirePolyLine(mergeType=IMPRINT, meshable= ON, points=((P1, P2), ))
		
		#Generate the name of the wire
		#str_wire = 'Wire-%d-Set-1' %(i-cnt_start+1)
		
		#Name the wire
		#mdb.models[modelName].parts[str_part].Set(edges=
    	#	mdb.models[modelName].parts[str_part].edges.getSequenceFromMask(('[#1 ]', ), ), name=str_wire)

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
        
#This function creates a CNT in Part module
def cnt_part(modelName, cnt_i, cnt_rad, cnt_start, cnt_end, cnt_coords):
    
    #Get the string for the CNT
    str_part = cnt_string_part(cnt_i)
    
    #Create a point to be able to generate the edges that will make the CNT
    mdb.models[modelName].Part(dimensionality=THREE_D, name=str_part, type=DEFORMABLE_BODY)
    mdb.models[modelName].parts[str_part].ReferencePoint(point=(0.0, 0.0, 0.0))
    
    #print("cnt_start=",cnt_start," cnt_end=",cnt_end)
    
    #Create edges, iterate over all points of the CNT
    generate_edges(modelName, cnt_start, cnt_end, cnt_coords, str_part)
    
    #Sweep an octagon along the edges
    Generate_Sweep(modelName, cnt_i, cnt_coords[cnt_end-1], cnt_coords[cnt_end], cnt_rad, cnt_start, cnt_end, str_part)
    
    #Delete the initial point as it is not used anymore
    del mdb.models[modelName].parts[str_part].features['RP']
	
#This function generates all CNT parts
def cnt_parts_all(modelName, N_CNTs, cnt_struct, cnt_coords):
	
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
		cnt_part(modelName, cnt_i, rad, acc_pts, acc_pts+N_p-1, cnt_coords)
        
		#Calculate the fraction of generated CNT parts
		frac = float(cnt_i)/float(N_CNTs)
		
		#Check if more than the fraction threshold has been generated
		if frac >= frac_thres:
		
			#Update the threshold
			while frac >= frac_thres:
			
				#Increase the threshold by a step
				frac_thres += check_step
			
			#Send message with percentage generated
			plog('   CNT parts generated: {:d} % ({} secs)\n'.format( int( (frac_thres - check_step)*100 ), time.time()-start_cnts))
		
		#Increase the number of accumulated points
		acc_pts += N_p

#This function creates the assembly by creating the instances of each part
def generate_assembly(modelName, N_CNTs, str_matrix):
	
	#Command necessary to create assembly
	mdb.models[modelName].rootAssembly.DatumCsysByDefault(CARTESIAN)
	
	#Generate name for matrix instance
	str_mat_inst = '%s-1' % (str_matrix)
	
	#Create instance for matrix
	mdb.models[modelName].rootAssembly.Instance(dependent=ON, name=str_mat_inst,
		part=mdb.models[modelName].parts[str_matrix])
	
	#Translate instance
	#i.e: endpoint - starting point = (P0[0], P0[1], P0[2]) - (P0[0]+matrixMeshSize, P0[1]+matrixMeshSize, matrixMeshSize)
	mdb.models[modelName].rootAssembly.translate(
	    instanceList=(str_mat_inst, ),
	    vector=(-matrixMeshSize, -matrixMeshSize, P0[2]-matrixMeshSize))
	
	#Create instances for each CNT
	for i in range(1, N_CNTs+1):
		
		#Generate name for CNT part
		str_cnt = cnt_string_part(i)
		
		#Generate name for CNT instance
		str_cnt_inst = cnt_string_instance(i)
		
		#Generate CNT instance
		mdb.models[modelName].rootAssembly.Instance(dependent=ON, name=str_cnt_inst,
			part=mdb.models[modelName].parts[str_cnt])		

#This function creates the matrix material and assigns it to a section
def materials_and_sections_matrix(modelName, str_matrix, matrixMaterial, matrixSection, matrixDensity, matrixModulus, matrixPoissonR):
	
	#Create matrix material
	mdb.models[modelName].Material(description='Polymer', name=matrixMaterial)
	mdb.models[modelName].materials[matrixMaterial].Elastic(table=((matrixModulus, matrixPoissonR), ))
	
	#Assign material to section
	mdb.models[modelName].HomogeneousSolidSection(material=matrixMaterial, name=matrixSection, thickness=None)
	
	#Assign section to matrix
	#print('cells=',len(mdb.models[modelName].parts[str_matrix].cells))
	mdb.models[modelName].parts[str_matrix].SectionAssignment(offset=0.0, offsetField='', offsetType=MIDDLE_SURFACE, 
		region=Region(cells=mdb.models[modelName].parts[str_matrix].cells), 
		sectionName=matrixSection, thicknessAssignment=FROM_SECTION)

#This function creates the CNT material and assigns it to a section
def materials_and_sections_cnt(modelName, N_CNTs, cntMaterial, cntSection, cntDensity, cntModulus, cntPoissonR):
	
	#Create CNT material
	mdb.models[modelName].Material(name=cntMaterial)
	mdb.models[modelName].materials[cntMaterial].Elastic(table=((cntModulus, cntPoissonR), ))
	
	#Assign material to section
	mdb.models[modelName].HomogeneousSolidSection(material=cntMaterial, name=cntSection, thickness=None)
	
	#Iterate over the CNTs
	for cnt_i in range(1, N_CNTs+1):
	
		#Get the string for the CNT part
		cnt_str = cnt_string_part(cnt_i)
		
		#Assign the CNT section to cnt_i
		mdb.models[modelName].parts[cnt_str].SectionAssignment(offset=0.0, offsetField='', offsetType=MIDDLE_SURFACE, 
			region=Region(cells=mdb.models[modelName].parts[cnt_str].cells.getSequenceFromMask(mask=('[#1 ]', ), )),
			sectionName=cntSection, thicknessAssignment=FROM_SECTION)

#This function creates an element set that contains the elements in the extended region of the RVE
#that need to be hidden in the visualization
def sets_for_elements_to_hide(modelName, matrixName, P0, Lxyz, matrixMeshSize, halfMatrixMeshSize):
	
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

#This function creates the sets that will be used when creating the embedded element constraints
def sets_for_embedded_elements(modelName, N_CNTs, str_matrix, str_host):
	
	#Set for the matrix
	mdb.models[modelName].rootAssembly.Set(
		cells=mdb.models[modelName].rootAssembly.instances[str_matrix + '-1'].cells,
		name=str_host)
	
	#Sets for the CNTs
	for cnt_i in range(1, N_CNTs+1):
	
		#Get the string for the CNT part
		cnt_str = cnt_string_part(cnt_i)
		
		#Generate name of set for cnt_i
		set_str = cnt_str + '_EESet'
		
		#Create the set for cnt_i
		mdb.models[modelName].rootAssembly.Set(
			cells=mdb.models[modelName].rootAssembly.instances[cnt_str + '-1'].cells,
			name=set_str)

#This functions creates the constraints for the embedded elements
def embedded_elements_constraints(modelName, N_CNTs, str_matrix, str_host):
	
	#Iterate over the CNTs
	for cnt_i in range(1, N_CNTs+1):
	
		#Get the string for the CNT part
		cnt_str = cnt_string_part(cnt_i)
		
		#Get the string for the CNT set
		set_str = cnt_str + '_EESet'
	
		#For cnt_i create an embedded element constraint
		mdb.models[modelName].EmbeddedRegion(
			absoluteTolerance=0.0, fractionalTolerance=0.05, toleranceMethod=BOTH, weightFactorTolerance=1e-06,
			embeddedRegion=mdb.models[modelName].rootAssembly.sets[set_str],
			hostRegion=mdb.models[modelName].rootAssembly.sets[str_host],
			name='EE-Constraint-%d' %(cnt_i))
	
#This function generates the mesh for the matrix
def generate_matrix_mesh(modelName, str_matrix, Lxyz, matrixMeshSize):
	
 # #Determine the size of the mesh for the polymer matrix
 # #First determine the size of the elements considering 10 elements per side
 # #In case the sample is not a cube, then first I need to find the shortest side
 # #This side will determine the size of the element
 # el_size_mat = min(Lxyz)/10.0
 #
 # #Check that the element size is at least twice the maximum CNT radius
 # if (el_size_mat - 2*cnt_rad_max <= Zero ):
 #
 # 	#The element size for the matrix is too small, so set it to be twice the maximum CNT radius
 # 	el_size_mat = 2*cnt_rad_max
 #

	#Mesh the matrix
	#deviationFactor and minSizeFactor have the default values from Abaqus
	mdb.models[modelName].parts[str_matrix].seedPart(deviationFactor=0.1, minSizeFactor=0.1, size=matrixMeshSize)
	mdb.models[modelName].parts[str_matrix].generateMesh()

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
        cnt_str = cnt_string_part(cnt_i)
        
        #Size of element
        cnt_el = 2.0*cnt_struct[cnt_i][1]

        #Mesh cnt_i, use its radius as the element size
        #deviationFactor and minSizeFactor have the default values from Abaqus
        mdb.models[modelName].parts[cnt_str].seedPart(
            deviationFactor=0.1, 
            minSizeFactor=0.1, 
            size=cnt_el)
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
            
            #Generete strings for messages
            str_try = 'dditional cut on part {} on point {}, cnt_start={} cnt_end={}\n'.format(str_part, j-cnt_start, cnt_start, cnt_end)
            str_except = 'Could not make a'+str_try
            str_try = 'A'+str_try
            
            #Cut the CNT cell at point j
            Partition_CNT_Cell_At_Point(modelName, cnt_rad, cnt_start, cnt_end, j, cnt_coords, str_part, hc, str_try, str_except)
    
    #Check if the invalid cuts can be done now
    if len(invalidCuts) > 0:  
        
        #Addtional cuts needed, iterate over those points
        for j in invalidCuts:
            
            #Generete strings for messages
            str_try = 'ut after second attempt on part {} on point {}, cnt_start={} cnt_end={}\n'.format(str_part, j-cnt_start, cnt_start, cnt_end)
            str_except = 'Could not make c'+str_try
            str_try = 'C'+str_try
            
            #Cut the CNT cell at point j
            Partition_CNT_Cell_At_Point(modelName, cnt_rad, cnt_start, cnt_end, j, cnt_coords, str_part, hc, str_try, str_except)

#This function partitions a cell at a given point
def Partition_CNT_Cell_At_Point(modelName, cnt_rad, cnt_start, cnt_end, idx, cnt_coords, str_part, hc, str_try, str_except):
            
    #Get the points of the CNT segments that share point i
    #End point of first CNT segment
    P1 = cnt_coords[idx-1]
    #Point shared by both CNT segments
    P2 = cnt_coords[idx]
    #End point of second CNT segment
    P3 = cnt_coords[idx+1]

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
        
        plog(str_try)
    except:
        
        plog(str_except)
        
#This function get the octagon needed to cut a CNT
def Get_octagon_to_cut_CNT(modelName, str_part, cnt_rad, P2, v1, v2, hc):
    
    #Get a unit vector that goes along the plane that bisects the angle 
    #between v1 and v2
    vm = make_unit((v1[0]+v2[0], v1[1]+v2[1], v1[2]+v2[2]))
    #print('vm')

    #Vector normal to v1, v2, and vm (all three lie on the same plane)
    N3 = cross(v2, v1)

    #Get the normal for the plane at the midpoint
    N = make_unit(cross(N3, vm))

    #Calculate a displacement to set the points of the cylinder
    disp = (N[0]*hc, N[1]*hc, N[2]*hc)

    #Calculate first point of Cylinder height
    C1 = (P2[0]+disp[0], P2[1]+disp[1], P2[2]+disp[2])

    #Calculate second point of Cylinder height
    C2 = (P2[0]-disp[0], P2[1]-disp[1], P2[2]-disp[2])

    #Calculate vector along the plane
    S = make_unit(cross(N3, v1))

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
    str_inst = cnt_string_instance(cnt_i)
    
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

#This functions creates a step and adds the boundary conditions for tension along the z-axis
def create_step_and_bcs(modelName, str_matrix, stpName, P0, corner, Lxyz):
	
	modelRoot = mdb.models[modelName].rootAssembly
	
	str_mat_instance = str_matrix + '-1'
    
	#Calculate half the lengths of the RVE along each direction
	half_x = Lxyz[0]*0.5
	half_y = Lxyz[1]*0.5
	half_z = Lxyz[2]*0.5
	
	#Set containing the 6 faces of the RVE (temperature will be applied to this set)
	#print('External_Faces')
	modelRoot.Set(
	    faces=modelRoot.instances[str_mat_instance].faces.findAt(
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
	    modelRoot.instances[str_mat_instance].faces.findAt(((corner[0], half_y, half_z), )), name='X_positive')
	modelRoot.Set(faces=
	    modelRoot.instances[str_mat_instance].faces.findAt(((P0[0], half_y, half_z), )), name='X_negative')
	modelRoot.Set(faces=
	    modelRoot.instances[str_mat_instance].faces.findAt(((half_x, corner[1], half_z), )), name='Y_positive')
	modelRoot.Set(faces=
	    modelRoot.instances[str_mat_instance].faces.findAt(((half_x, P0[1], half_z), )), name='Y_negative')
	modelRoot.Set(faces=
	    modelRoot.instances[str_mat_instance].faces.findAt(((half_x, half_y, corner[2]), )), name='Z_positive')
	modelRoot.Set(faces=
	    modelRoot.instances[str_mat_instance].faces.findAt(((half_x, half_y, P0[2]), )), name='Z_negative')
	
	#############################################################################
	#Equations with the same reference point (RP_#) terms will be added
	#E.g.: 1*X_positive + 0.5*RP_X + (-1*X_negative) + 0.5*RP_X = X_Positive - X_Negative + RP_X = 0
	#If RPs are needed, then this sould be enabled
	
	modelRoot = mdb.models[modelName]
	
	modelRoot.Equation(name='EQ_X_Positive', terms=((1.0, 'X_positive', 1), (-0.5, 'RP_X', 1)))
	modelRoot.Equation(name='EQ_X_Negative', terms=((-1.0, 'X_negative', 1), (-0.5, 'RP_X', 1)))
	
	modelRoot.Equation(name='EQ_Y_Positive', terms=((1.0, 'Y_positive', 2), (-0.5, 'RP_Y', 2)))
	modelRoot.Equation(name='EQ_Y_Negative', terms=((-1.0, 'Y_negative', 2), (-0.5, 'RP_Y', 2)))
	
	modelRoot.Equation(name='EQ_Z_Positive', terms=((1.0, 'Z_positive', 3), (-0.5, 'RP_Z', 3)))
	modelRoot.Equation(name='EQ_Z_Negative', terms=((-1.0, 'Z_negative', 3), (-0.5, 'RP_Z', 3)))
	
	#############################################################################
	#Set name of the amplitude
	ampName = 'Amp-1'
	
	#Create a simple amplitude
	mdb.models[modelName].TabularAmplitude(
	    data=((0.0, 0.0), (1.0, 1.0)), 
	    name=ampName, 
	    smooth=SOLVER_DEFAULT, 
	    timeSpan=STEP)
	
	#############################################################################
	#Create mechanical step
	mdb.models[modelName].StaticStep(
	    initialInc=inicialIncrement, 
	    name=stpName, 
	    noStop=OFF, 
	    previous='Initial', 
	    timeIncrementationMethod=FIXED)
	
	#############################################################################
	#Name for the BC
	bcDispName = 'BC-Disp'
	        
	#Generate an array for the reference points to which PBCs will be applied
	pointsPBCs = []
	
	#Append indicated by the flags
	#Check if displacement along the x-direction is applied
	if dispXflag:
	    
	    #Append reference point to apply displacement along the x-direction
	    pointsPBCs.append(mdb.models[modelName].rootAssembly.sets['RP_X'].referencePoints[0])
	
	#Check if displacement along the y-direction is applied
	if dispYflag:
	    
	    #Append reference point to apply displacement along the y-direction
	    pointsPBCs.append(mdb.models[modelName].rootAssembly.sets['RP_Y'].referencePoints[0])
	
	#Check if displacement along the z-direction is applied
	if dispZflag:
	    
	    #Append reference point to apply displacement along the z-direction
	    pointsPBCs.append(mdb.models[modelName].rootAssembly.sets['RP_Z'].referencePoints[0])
	
	#Apply the displacement BC to 
	mdb.models[modelName].DisplacementBC(
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
	    mdb.models[modelName].boundaryConditions[bcDispName].setValues(u1=dispX)
	
	#Check if displacement along the y-direction is applied
	if dispYflag:
	    
	    #Se the displacement BC along the y-direction
	    mdb.models[modelName].boundaryConditions[bcDispName].setValues(u2=dispY)
	
	#Check if displacement along the z-direction is applied
	if dispZflag:
	    
	    #Se the displacement BC along the z-direction
	    mdb.models[modelName].boundaryConditions[bcDispName].setValues(u3=dispZ)
	    
	#############################################################################
	#Set the output request
	mdb.models[modelName].fieldOutputRequests['F-Output-1'].setValues(variables=('S', 'E', 'U'))
	
#This function creates two sets for two of the vertices of the matrix (sample), where each set has only one node
#These nodes are on the diagonal of the cuboid that defines the matrix (sample) 
def create_sets_for_matrix(modelName, P0, Lxyz, str_matrix):
	
	#Create the set of the lower left corner
	#NOTE: The parenthesis, comma and the space in the nodes argument is beacuse a tuple is needed but
	#the operation getClosest returns an object. The parenthesis, comma and space make it a tuple
	mdb.models[modelName].rootAssembly.Set(
		nodes=mdb.models[modelName].rootAssembly.instances[str_matrix+"-1"].nodes.getByBoundingSphere(center=P0, radius=0.001),
		name='Matrix0')
	
	#print('P0=',P0,' str_matrix=',str_matrix)
	#nd = mdb.models[modelName].rootAssembly.instances[str_matrix+"-1"].nodes.getByBoundingSphere(center=P0, radius=0.001)
	#print('nd0=',nd[0].coordinates)
	
	#Create opposite corner
	corner = (P0[0]+Lxyz[0], P0[1]+Lxyz[1], P0[2]+Lxyz[2])
	
	#Create the set of the opposite corner
	mdb.models[modelName].rootAssembly.Set(
		nodes=mdb.models[modelName].rootAssembly.instances[str_matrix+"-1"].nodes.getByBoundingSphere(center=corner, radius=0.001),
		name='Matrix1')
	#nd = mdb.models[modelName].rootAssembly.instances[str_matrix+"-1"].nodes.getByBoundingSphere(center=corner, radius=0.001)
	#print('nd1=',nd[0].coordinates)
		
	#Perform union of two sets
	#mdb.models[modelName].rootAssembly.SetByBoolean(
	#	name='Matrix2', 
	#	operation=UNION, 
	#	sets=(mdb.models[modelName].rootAssembly.sets['Matrix0'], mdb.models[modelName].rootAssembly.sets['Matrix1'])
	#	)

#This function creates all sets for the nodes that correspond to the centerline of that CNT
def Create_All_Sets_For_CNT_Points(modelName, N_CNTs, cnt_struct, cnt_coords):

    #Initializae the number of accumulated points to zero
    acc_pts = 0

    #Iterate over all CNTs
    for cnt_i in range(1, N_CNTs+1):

        #Number of points in CNTi and its radius
        N_p = int(cnt_struct[cnt_i][0])
        cnt_rad = cnt_struct[cnt_i][1]

        #Get the strings for part and instance corresponding to cnt_i
        cnt_str_part = string_part('CNT', cnt_i)
        cnt_str_inst = string_instance('CNT', cnt_i)
        
        #Get CNT start and end
        cnt_start = acc_pts
        cnt_end = acc_pts+N_p-1
        
        #Make set with the CNT points
        Create_Set_For_CNT_Points(modelName, cnt_i, cnt_rad, cnt_start, cnt_end, cnt_coords, cnt_str_inst)

        #Get the name of the node set
        node_set_str = cnt_string_node_set(cnt_i)
            
        #Get the number of nodes in the set
        n_nodes = len(mdb.models[modelName].rootAssembly.sets[node_set_str].nodes)
        #plog('n_nodes = {}\n'.format(n_nodes))
        
        #Check that the set has the right amount if nodes
        if n_nodes != N_p:
            
            if n_nodes == 0:
                plog('Empty set due to irregular mesh: {}\n'.format(node_set_str))
            else:
                plog('Set {} does not have the expected number of nodes ({} nodes but {} points) likely due to an irregular mesh\n'.format(node_set_str,n_nodes,N_p))
            
            #Make two extra cuts
            plog('Attempting two more cuts on CNT cell...\n')
            
            #Delete intial set
            del mdb.models[modelName].rootAssembly.sets[node_set_str]
            
            #Calculate height of cylinder
            hc = cnt_rad*0.1
            
            #Generate strings for messages
            str_try = 'dditional cut (to generate set) on part {} on point {}, cnt_start={} cnt_end={}\n'.format(cnt_str_part, 1, cnt_start, cnt_end)
            str_except = 'Could not make a'+str_try
            str_try = 'A'+str_try
            
            #Cut the CNT cell at second point in the CNT (cnt_start+1)
            Partition_CNT_Cell_At_Point(modelName, cnt_rad, cnt_start, cnt_end, cnt_start+1, cnt_coords, cnt_str_part, hc, str_try, str_except)
            
            #Generate strings for messages
            str_try = 'dditional cut (to generate set) on part {} on point {}, cnt_start={} cnt_end={}\n'.format(cnt_str_part, N_p-2, cnt_start, cnt_end)
            str_except = 'Could not make a'+str_try
            str_try = 'A'+str_try
            
            #Cut the CNT cell at previous to last point in CNT (cnt_end-2)
            Partition_CNT_Cell_At_Point(modelName, cnt_rad, cnt_start, cnt_end, cnt_end-2, cnt_coords, cnt_str_part, hc, str_try, str_except)
            
            #Mesh again
            mdb.models[modelName].parts[cnt_str_part].seedPart(
                deviationFactor=0.1, 
                minSizeFactor=0.1, 
                size=2.0*cnt_rad)
            mdb.models[modelName].parts[cnt_str_part].generateMesh()
        
            #This seems to be required by Abaqus to update the mesh
            mdb.models[modelName].rootAssembly.regenerate()
            
            #Create set again
            Create_Set_For_CNT_Points(modelName, cnt_i, cnt_rad, cnt_start, cnt_end, cnt_coords, cnt_str_inst)
            
            #Get the number of nodes in the set
            n_nodes = len(mdb.models[modelName].rootAssembly.sets[node_set_str].nodes)
            #plog('updated n_nodes = {}\n'.format(n_nodes))

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
    cntSet = mdb.models[modelName].rootAssembly.Set(
        nodes=allNodes,
        name=node_set_str)

    #Print the length of the set
    #plog('%s nodes=%d points=%d\n'%(node_set_str, len(mdb.models[modelName].rootAssembly.sets[node_set_str].nodes), cnt_end+1-cnt_start))
    
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
        #plog('%s nodes=%d points=%d\n'%(node_set_str, len(mdb.models[modelName].rootAssembly.sets[node_set_str].nodes), cnt_end+1-cnt_start))
	
	
######################################---ABAQUS FUNCTIONS---########################################
######################################----------END---------########################################

######################################---MAIN PROGRAM---########################################

#Time for whole Abaqus model
start0 = time.time()

#Read the csv file and turn it into a list
#By using the flag QUOTE_NONNUMERIC, the reader assumes there are
#floats unless the content in the csv file has single quotes
#Thus, conversion from string to float is made automatically

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

#Read the sample geometry
with open(sample_file) as f:
    (P0, Lxyz) = list(csv.reader(f, quoting=csv.QUOTE_NONNUMERIC))
#print(P0)
#print(Lxyz)

#Determine the size of the mesh for the polymer matrix
#As a default mesh size, consider 10 elements per side
#But elment size should be at least twide the radius
matrixMeshSize = max(min(Lxyz)/elementsPerSide, 2.5*cnt_rad_max)
#print('matrixMeshSize=',matrixMeshSize)
halfMatrixMeshSize = 0.5*matrixMeshSize

#Calculate furthest corner of RVE
corner = (P0[0]+Lxyz[0], P0[1]+Lxyz[1], P0[2]+Lxyz[2])
#print(corner)

#Calculate lengths of extended RVE
Lxyz_ext = (Lxyz[0] + 2.0*matrixMeshSize, Lxyz[1] + 2.0*matrixMeshSize, Lxyz[2] + 2.0*matrixMeshSize)

#Generate the matrix part
matrix_part(modelName, P0, Lxyz, Lxyz_ext, matrixName)
#mdb.models[modelName].parts[matrixName].DatumPointByCoordinate(corner)
    
#Get the number of CNTs
N_CNTs = int(cnt_struct[0][0])
#A small fixed number used for debugging and testing
#N_CNTs=2

#Name of the job to be used based on its parameters
#CNT-'Number of CNTs in the RVE'
#EPS-'Number of elements per side'
jobName = 'CNT-'+str(N_CNTs)+'_EPS-'+str(elementsPerSide)

#Name of the file to save print messages
print_file = jobName + '.txt'
pfile = open(print_file, "a")

plog('There are ' + str(N_CNTs) + ' CNTs inside the RVE.\n')

start = time.time()

#Generate all CNT parts
cnt_parts_all(modelName, N_CNTs, cnt_struct, cnt_coords)

end = time.time()
plog("Time for part generation: {}\n".format(end-start))

#Generate materials and assign sections
materials_and_sections_matrix(modelName, matrixName, matrixMaterial, matrixSection, matrixDensity, matrixModulus, matrixPoissonR)
materials_and_sections_cnt(modelName, N_CNTs, cntMaterial, cntSection, cntDensity, cntModulus, cntPoissonR)

start = time.time()
#Generate assembly
generate_assembly(modelName, N_CNTs, matrixName)

end = time.time()
plog("Time for instance generation: {}\n".format(end-start))

#Create sets that will be used when creating the embedded element constraints
sets_for_embedded_elements(modelName, N_CNTs, matrixName, str_host)

#Create embedded elements constraints
embedded_elements_constraints(modelName, N_CNTs, matrixName, str_host)

#Create Step and add boundary conditions
create_step_and_bcs(modelName, matrixName, stpName, P0, corner, Lxyz)

start = time.time()
#Generate meshes
generate_matrix_mesh(modelName, matrixName, Lxyz, matrixMeshSize)
Generate_CNT_Meshes(modelName, N_CNTs, cnt_struct, cnt_coords)

end = time.time()
plog("Time for meshing: {}\n".format(end-start))
start = time.time()
#Create set for elements to hide in visualization
sets_for_elements_to_hide(modelName, matrixName, P0, Lxyz, matrixMeshSize, halfMatrixMeshSize)

#Create sets for the matrix (sample in the C++ code)
#NOTE: Sets are generated on root assembly
create_sets_for_matrix(modelName, P0, Lxyz, matrixName)

#Create sets for central CNT nodes
#NOTE: Sets are generated on root assembly
Create_All_Sets_For_CNT_Points(modelName, N_CNTs, cnt_struct, cnt_coords)

end = time.time()
plog("Time for creating CNT point sets and elemnts to hide: {}\n".format(end-start))

if createJob == 1:
    
    #Create and submit job using Abaqus default values
    mdb.Job(atTime=None, contactPrint=OFF, description='', echoPrint=OFF, 
    	explicitPrecision=SINGLE, getMemoryFromAnalysis=True, historyPrint=OFF, 
    	memory=90, memoryUnits=PERCENTAGE, model=modelName, modelPrint=OFF, 
    	multiprocessingMode=DEFAULT, name=jobName, nodalOutputPrecision=SINGLE, 
    	numCpus=1, numGPUs=0, queue=None, resultsFormat=ODB, scratch='', type=
    	ANALYSIS, userSubroutine='', waitHours=0, waitMinutes=0)
    
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
        plog("Data check on input file: {} secs.\n".format(time.time()-start))

    if submitJob == 1:

        #Submit job and save time of submission
        start = time.time()
        mdb.jobs[jobName].submit(consistencyChecking=OFF)
        
        #Make the Python script to wait for the Abaqus job to finish
        #In this way the script can measure the execution time of the Abaqus simulation
        mdb.jobs[jobName].waitForCompletion()

        #Calcualte job executionand simulation time
        plog("Time for Job execution: {}\n".format(time.time()-start))
        
plog("Time for Abaqus model: {}\n".format(time.time()-start0))
pfile.close()
