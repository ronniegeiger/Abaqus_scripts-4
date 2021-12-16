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

######################################---GLOBAL VARIABLES---########################################

#Files to read
coord_file = 'C:\Users\grafenicas\Documents\CNTs_embedde_elements\Job-testing-mix\cnt_coordinates.csv'
struct_file = 'C:\Users\grafenicas\Documents\CNTs_embedde_elements\Job-testing-mix\cnt_struct.csv'
sample_file = 'C:\Users\grafenicas\Documents\CNTs_embedde_elements\Job-testing-mix\sample_geom.csv'
gs_file = 'C:\Users\grafenicas\Documents\CNTs_embedde_elements\Job-testing-mix\gnp_data.csv'
	
#Define the string for the name of the matrix part
str_matrix = 'Matrix'

#String for the host set (i.e., the matrix)
str_host = 'host_Set'

#Define the sheetsize for the drawing space
sheetSize = 20.0
#Mesh size for the matrix (um)
matrixMeshSize = 0.1
#Define the EmbeddedMesh/HostMesh ratio
meshRatio = 0.5
#Temperature difference (C)
tempApplied = 50.0

#Cosine of 45 degrees (PI/4)
cos45 = 0.7071067812

#Factor to convert from radians to degrees
rad_to_deg = 180/math.pi

#"Zero" for comparing floating point numbers
Zero = 1e-7

#Displacement to be applied (microns)
disp = -0.05

#Matrix properties
#Define the name and section of the matrix
matrixMaterial = 'Matrix_mat'
matrixSection = 'Matrix_sec'
#Define the mass density
matrixDensity = 1235
#Define the elastic modulus (MPa)
matrixModulus = 2500
#Define the Poisson ratio
matrixPoissonR = 0.37
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
cntMaterial = 'CNT_mat'
cntSection = 'CNT_sec'
#Define the mass density (kg/m3)
cntDensity = 1800
#Define the elastic modulus (GPa)
cntModulus = 1e6
#Define the Poisson ratio
cntPoissonR = 0.35
#Define de coefficient of thermal expansion (e-5 1/C)
cntExpCoeff = 0.5
#Define the thermal conductivity (W/m*K)
cntThermConductivity = 3000
#Define the specific heat (J/mol*K)
cntSpecHeat = 7
#Define the electrical conductivity (S/m)
cntElecConductivity = 1e7
#Maximum CNT radius
cnt_rad_max = 0.005


#Define the name of the filler
gsMaterial = 'Graphene'
gsSection = 'Graphene-sec'
#Define the mass density (kg/m3)
gsDensity = 2000
#Define the elastic modulus (MPa)
gsModulus = 1e6
#Define the Poisson ratio
gsPoissonR = 0.165
#Define de coefficient of thermal expansion (e-5 1/C)
gsExpCoeff = 0.5
#Define the thermal conductivity (W/m*K)
gsThermConductivity = 3000
#Define the specific heat (J/mol*K)
gsSpecHeat = 7
#Define the electrical conductivity (S/m)
gsElecConductivity = 1e7

#Cutting box name
biggerBoxName = 'BIGGER_BOX'
#Name of the analysis model (unchangeable)
modelName = 'Model-1'

#Analsys type for the element
# 1 = coupled temperature-displacement
# 2 = coupled thermal-electrical-structural
# 3 = thermal-electric
# 4 = only displacement
elemTypeCode = 4

#Different element types
elemTemp_Disp = [C3D8T, C3D6T, C3D4T]
elemTemp_Elec_Disp = [Q3D8, Q3D6, Q3D4]
elemTemp_Elec = [DC3D8E, DC3D6E, DC3D4E]
elem_Disp = [C3D8R, C3D6, C3D4]
selectedElementCode = []

#Selecting a desired element type
if elemTypeCode == 1:
    selectedElementCode = elemTemp_Disp
elif elemTypeCode == 2:
    selectedElementCode = elemTemp_Elec_Disp
elif elemTypeCode == 3:
    selectedElementCode = elemTemp_Elec
elif elemTypeCode == 4:
    selectedElementCode = elem_Disp
else:
    print('WARNING. Invalid element type. Elment type code was set to displacement elements (elemTemp_Disp)')
    selectedElementCode = elem_Disp

#Embedded element mesh size
eeMeshSize = matrixMeshSize*meshRatio

######################################---GLOBAL VARIABLES---########################################

######################################---ABAQUS FUNCTIONS---########################################

#This function generates a string for the filler part
def string_part(filler, i):
    #Return the string in the format FILLER-i
    return '%s-%d' %(filler, i)

#This function generates a string for the filler instance
def string_instance(filler, i):
    #Return the string in the format FILLER-i-1
    return '%s-%d-1' %(filler, i)
    
def ee_string(filler, i):
    return 'EE-%s-%d'%(filler, i)
   
#This function generates a string for the node set of filler i
def string_node_set(filler, i):
	#Return the string in the format CNT-nodes-cnt_i
	return '%s-%d-nodes' %(filler, i)
   
def matrix_part(P0, Lxyz, str_matrix, sheetSize):
	
	#Create a sketch
	mdb.models['Model-1'].ConstrainedSketch(name='__profile__', sheetSize=sheetSize)
	
	#Define the two corners of a rectangle on the xy plane
	mdb.models['Model-1'].sketches['__profile__'].rectangle(
		point1=(P0[0], P0[1]), point2=(P0[0]+Lxyz[0], P0[1]+Lxyz[1]))
		
	#Name the part
	mdb.models['Model-1'].Part(dimensionality=THREE_D, name=str_matrix, type= DEFORMABLE_BODY)
	
	#Use the length along the z-direction as the extrusion depth
	mdb.models['Model-1'].parts[str_matrix].BaseSolidExtrude(depth=Lxyz[2], sketch=
	    mdb.models['Model-1'].sketches['__profile__'])
	    
	#Delete the sketch
	del mdb.models['Model-1'].sketches['__profile__']
	
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
	
#This function generates all edges of a CNT
def generate_edges(cnt_start, cnt_end, cnt_coords, str_part):
	
	#Create edges, iterate over all points of the CNT
	for i in range(cnt_start, cnt_end):
	#for i in range(cnt_start, cnt_start+11):
		
		#First point of the CNT segment
		P1 = cnt_coords[i]
		#Second point of the CNT segment
		P2 = cnt_coords[i+1]
		
		#Generate the segment
		mdb.models['Model-1'].parts[str_part].WirePolyLine(mergeType=IMPRINT, meshable= ON, points=((P1, P2), ))
		
		#Generate the name of the wire
		#str_wire = 'Wire-%d-Set-1' %(i-cnt_start+1)
		
		#Name the wire
		#mdb.models['Model-1'].parts[str_part].Set(edges=
    	#	mdb.models['Model-1'].parts[str_part].edges.getSequenceFromMask(('[#1 ]', ), ), name=str_wire)

#This function generates the sweeped part
def generate_sweep(P1, P2, cnt_rad, cnt_start, cnt_end, str_part):
	#The sketching plane is perpendicular to the last segment
	#The last segment has the points P1 and P2, so the z-axis of the sketch plane is 
	#aligned to this last segment and goes in the direction from P1 to P2
	#A transformation matrix is generated to align the z-axis to the last segment
	R = rotation_matrix(P1, P2)
	
	#Create the sketching plane using the rotation matrix and the last point in the CNT
	mdb.models['Model-1'].ConstrainedSketch(gridSpacing=0.001, name='__profile__', sheetSize=0.076, transform=(
	    R[0][0], R[0][1], R[0][2],
	    R[1][0], R[1][1], R[1][2], 
	    R[2][0], 0.0, R[2][2], 
	    P2[0], P2[1], P2[2]))
	mdb.models['Model-1'].sketches['__profile__'].sketchOptions.setValues( decimalPlaces=5)
	mdb.models['Model-1'].sketches['__profile__'].ConstructionLine(point1=(-0.038, 0.0), point2=(0.038, 0.0))
	mdb.models['Model-1'].sketches['__profile__'].ConstructionLine(point1=(0.0,-0.038), point2=(0.0, 0.038))
	mdb.models['Model-1'].parts[str_part].projectReferencesOntoSketch(
		filter=COPLANAR_EDGES, sketch=mdb.models['Model-1'].sketches['__profile__'])
	
	#Calculate the radius multiplied by cos(PI/4)
	new_rad = cnt_rad*cos45
	
    #Construct a regular octagon
    #Vertex 1-2
	mdb.models['Model-1'].sketches['__profile__'].Line(point1=(0.0, cnt_rad), point2=(new_rad, new_rad))
	#Vertex 2-3
	mdb.models['Model-1'].sketches['__profile__'].Line(point1=(new_rad, new_rad), point2=(cnt_rad, 0.0))
	#Vertex 3-4
	mdb.models['Model-1'].sketches['__profile__'].Line(point1=(cnt_rad, 0.0), point2=(new_rad, -new_rad))
	#Vertex 4-5
	mdb.models['Model-1'].sketches['__profile__'].Line(point1=(new_rad,-new_rad), point2=(0.0, -cnt_rad))
	#Vertex 5-6
	mdb.models['Model-1'].sketches['__profile__'].Line(point1=(0.0, -cnt_rad), point2=(-new_rad, -new_rad))
	#Vertex 6-7
	mdb.models['Model-1'].sketches['__profile__'].Line(point1=(-new_rad, -new_rad), point2=(-cnt_rad, 0.0))
	#Vertex 7-8
	mdb.models['Model-1'].sketches['__profile__'].Line(point1=(-cnt_rad, 0.0), point2=(-new_rad, new_rad))
	#Vertex 8-1
	mdb.models['Model-1'].sketches['__profile__'].Line(point1=(-new_rad, new_rad), point2=(0.0, cnt_rad))	
	
	#Select the last edge, which has number equal to (number edges-1)
	#The number of edges is equal to (number of points-1)
	#The number of points is (cnt_end-cnt_start+1)
	#Then, the las edge has number ((cnt_end-cnt_start+1)-1-1)=(cnt_end-cnt_start-1)
	mdb.models['Model-1'].parts[str_part].SolidSweep(
		path=mdb.models['Model-1'].parts[str_part].edges, 
		profile=mdb.models['Model-1'].sketches['__profile__'], 
		sketchOrientation=RIGHT, 
		sketchUpEdge=mdb.models['Model-1'].parts[str_part].edges[cnt_end-cnt_start-1])
	    
	#Delete sketch
	del mdb.models['Model-1'].sketches['__profile__']
	
	
#This function retunrs a unit vector going from P1 towards P2
def get_unit_vector(P1, P2):
	
	#Get the vector going from P1 towards P2
	v = (P2[0]-P1[0], P2[1]-P1[1], P2[2]-P1[2])
	
	#Calculate the length of vector v
	length = sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2])
	
	#Return v as a unit vector
	return (v[0]/length, v[1]/length, v[2]/length)
	
#Cross product of two vectors
def cross(v1, v2):
	
	#Calculate the compoenents of the cross product
	x = v1[1]*v2[2] - v1[2]*v2[1]
	y = -v1[0]*v2[2] + v1[2]*v2[0]
	z = v1[0]*v2[1] - v1[1]*v2[0]
	
	#Return the cross product
	return (x,y,z)
	
#This function generates all CNT parts
def cnt_parts_all(N_CNTs, cnt_struct, cnt_coords):
	
	#Number of accumulated points
	acc_pts = 0
	
	#Iterate over the number of CNTs
	for cnt_i in range(1, N_CNTs+1):
		
		#Number of points in CNTi and its radius
		N_p = int(cnt_struct[cnt_i][0])
		rad = cnt_struct[cnt_i][1]
		
		#Create all CNTs, one part per CNT
		#The first point is given by the accumulated number of points
		#The last point is the accumulated number of points plus the number of points of CNTi minus 1
		cnt_part(cnt_i, rad, acc_pts, acc_pts+N_p-1, cnt_coords)
		
		#Increase the number of accumulated points
		acc_pts += N_p
    		    
#This function creates a CNT in Part module
def cnt_part(cnt_i, cnt_rad, cnt_start, cnt_end, cnt_coords):
	
	#Get the string for the CNT
	str_part = string_part('CNT', cnt_i)
	
	#Create a point to be able to generate the edges that will make the CNT
	mdb.models['Model-1'].Part(dimensionality=THREE_D, name=str_part, type=DEFORMABLE_BODY)
	mdb.models['Model-1'].parts[str_part].ReferencePoint(point=(0.0, 0.0, 0.0))
	
	#print("cnt_start=",cnt_start," cnt_end=",cnt_end)
	
	#Create edges, iterate over all points of the CNT
	generate_edges(cnt_start, cnt_end, cnt_coords, str_part)
	
	#Sweep an octagon along the edges
	generate_sweep(cnt_coords[cnt_end-1], cnt_coords[cnt_end], cnt_rad, cnt_start, cnt_end, str_part)
	
	#Delete the initial point as it is not used anymore
	del mdb.models['Model-1'].parts[str_part].features['RP']
		
	#datum_points(cnt_rad, cnt_start, cnt_end, cnt_coords, str_part)

#Create graphene sheet
def gs_parts_all(N_GSs, gs_data, modelName, sheetSize):
    
    #Loop to define the number of graphene sheets to create
    for i in range(N_GSs):
        
        #Create a graphene sheet
        Create_GS(gs_data[i], modelName, sheetSize, i)
        #Create_NodeSet_Gs(modelName, numPart)    
    
def Create_GS(GS_geometry, model, sheetsz, numPart):
    
    #Get the length along each direction
    length_x = GS_geometry[0]
    length_y = GS_geometry[1]
    thickness = GS_geometry[2]
    
    #Calculate half length along x and y directions
    half_x = length_x/2.0
    half_y = length_y/2.0
    
    #Get the part name
    str = string_part('GS', numPart)

    mdb.models[model].ConstrainedSketch(name='__profile__', sheetSize=sheetsz)
    mdb.models[model].sketches['__profile__'].rectangle(
        point1=(-half_x, -half_y), 
        point2=( half_x,  half_y))
    mdb.models[model].Part(dimensionality=THREE_D, name=str, type=DEFORMABLE_BODY)
    mdb.models[model].parts[str].BaseSolidExtrude(depth=thickness, sketch=mdb.models[model].sketches['__profile__'])
    del mdb.models[model].sketches['__profile__']
    
#This function creates the assembly by creating the instances of each part
def generate_assembly(N_CNTs, N_GSs, model, str_matrix, gs_data):
	
	#Necessary command to create assembly
	mdb.models[model].rootAssembly.DatumCsysByDefault(CARTESIAN)
	
	#Generate name for matrix instance
	str_mat_inst = '%s-1' % (str_matrix)
	
	#Create instance for matrix
	mdb.models[model].rootAssembly.Instance(dependent=ON, name=str_mat_inst,
		part=mdb.models[model].parts[str_matrix])
		
	if N_CNTs != 0:
		#Create the instances for all CNTs
		generate_cnt_instances(model, N_CNTs)
	
	if N_GSs != 0:	
		#Create the cutting box using one particular GS as reference
		#Create_CuttingBox(modelName, str_matrix, biggerBoxName, gs_data[1])
		
		#Create the instances for all GSs
		generate_gs_instances(model, N_GSs, gs_data)
	
def generate_cnt_instances(model, N_CNTs):
	
	#Create instances for each CNT
	for i in range(1, N_CNTs+1):
		
		#Generate name for CNT part
		str_cnt = string_part('CNT', i)
		
		#Generate name for CNT instance
		str_cnt_inst = string_instance('CNT', i)
		
		#Generate CNT instance
		mdb.models[model].rootAssembly.Instance(dependent=ON, name=str_cnt_inst,
			part=mdb.models[model].parts[str_cnt])

def generate_gs_instances(model, N_GSs, gs_data):
    
    #Create instances for each GS
    for i in range(N_GSs):
        
        #Generate name for part
        str_part = string_part('GS', i)
        
        #Generate name for instance
        str_inst = string_instance('GS', i)
        
        #Generate part instance
        mdb.models[model].rootAssembly.Instance(dependent=ON, name=str_inst,
            part=mdb.models[model].parts[str_part])
        
        #Apply translation and rotation to GS to get its final position 
        Translate_and_Rotate_GS(model, i, str_inst, gs_data[i])
        
        #Check if GS is partially outside the sample and cut it if so
        #Cut_GS(model, i, str_inst, str_part)
	
#Rotate the graphene sheets
def Translate_and_Rotate_GS(model, i, str_inst, GS_parameters):
	
	#Get the paramters for translation and rotation (for readability of the script)
	thickness = GS_parameters[2]
	y_angle = GS_parameters[3]*rad_to_deg
	z_angle = GS_parameters[4]*rad_to_deg
	final_point = GS_parameters[5:]

	#Translate GS centroid (x0,y0,z0) to the origin
	#i.e: endpoint - starting point = (0,0,0) - (x0,y0,z0)
	mdb.models[model].rootAssembly.translate(
		instanceList=(str_inst, ), 
		vector=(0.0, 0.0, -thickness/2))

	#Rotation in y-axis
	mdb.models[model].rootAssembly.rotate(
		angle = y_angle, 
		axisDirection=(0.0, 1.0, 0.0), 
		axisPoint=(0.0, 0.0, 0.0), 
		instanceList=(str_inst, ))

	#Rotation in z-axis
	mdb.models[model].rootAssembly.rotate(
		angle = z_angle, 
		axisDirection=(0.0, 0.0, 1.0), 
		axisPoint=(0.0, 0.0, 0.0), 
		instanceList=(str_inst, ))

	#Translate GS centroid (currently in the origin) to its new position (x1,y1,z1)
	#I.e: endpoint - starting point = (x1,y1,z1) - (0,0,0)
	mdb.models[model].rootAssembly.translate(
		instanceList=(str_inst, ), 
		vector=final_point)

	# mdb.models['Model-1'].parts['GS-0'].Set(name='SET-GS-0', nodes=mdb.models['Model-1'].parts['GS-0'].nodes.getSequenceFromMask(mask=(
	#	 '[#20301 #0 #201000 #0 #3010000 #2 ]', ), ))

#This function creates the matrix material and assigns it to a section
def materials_and_sections_matrix(model, str_matrix, matrixMaterial, matrixSection, matrixDensity, matrixModulus, matrixPoissonR):
	
	#Create matrix material
	mdb.models[model].Material(description='Polymer', name=matrixMaterial)
	mdb.models[model].materials[matrixMaterial].Elastic(table=((matrixModulus, matrixPoissonR), ))
	
	#Assign material to section
	mdb.models[model].HomogeneousSolidSection(material=matrixMaterial, name=matrixSection, thickness=None)
	
	#Assign section to matrix
	#print('cells=',len(mdb.models['Model-1'].parts[str_matrix].cells))
	mdb.models[model].parts[str_matrix].SectionAssignment(offset=0.0, offsetField='', offsetType=MIDDLE_SURFACE, 
		region=Region(cells=mdb.models[model].parts[str_matrix].cells), 
		sectionName=matrixSection, thicknessAssignment=FROM_SECTION)

#This function creates the CNT material and assigns it to a section
def materials_and_sections_cnt(model, N_CNTs, cntMaterial, cntSection, cntDensity, cntModulus, cntPoissonR):
	
	#Create CNT material
	mdb.models[model].Material(name=cntMaterial)
	mdb.models[model].materials[cntMaterial].Elastic(table=((cntModulus, cntPoissonR), ))
	
	#Assign material to section
	mdb.models[model].HomogeneousSolidSection(
		material=cntMaterial, 
		name=cntSection, 
		thickness=None)
	
	#Iterate over the CNTs
	for cnt_i in range(1, N_CNTs+1):
	
		#Get the string for the CNT part
		cnt_str = string_part('CNT', cnt_i)
		
		#Assign the CNT section to cnt_i
		mdb.models[model].parts[cnt_str].SectionAssignment(offset=0.0, offsetField='', offsetType=MIDDLE_SURFACE, 
			region=Region(cells=mdb.models[model].parts[cnt_str].cells.getSequenceFromMask(mask=('[#1 ]', ), )),
			sectionName=cntSection, thicknessAssignment=FROM_SECTION)

#This function creates the GS material and assigns it to a section
def materials_and_sections_gs(model, N_GSs, gsMaterial, gsSection, gsDensity, gsModulus, gsPoissonR):
    
    #Create GS material
    mdb.models[model].Material(name=gsMaterial)
    mdb.models[model].materials[gsMaterial].Elastic(table=((gsModulus, gsPoissonR), ))
    
    #Assign material to section
    mdb.models[model].HomogeneousSolidSection(
        material=gsMaterial, 
        name=gsSection, 
        thickness=None)
    
    #Iterate over the CNTs
    for gs_i in range(N_GSs):
    
        #Get the string for the GS part
        gs_str = string_part('GS', gs_i)
        
        #Assign the CNT section to gs_i
        mdb.models[model].parts[gs_str].SectionAssignment(
            offset=0.0, offsetField='', offsetType=MIDDLE_SURFACE, 
            region=Region(cells=mdb.models[model].parts[gs_str].cells.getSequenceFromMask(mask=('[#1 ]', ), )),
            sectionName=gsSection, thicknessAssignment=FROM_SECTION)
    
#This function creates the sets that will be used when creating the embedded element constraints
def sets_for_embedded_elements(P0, Lxyz, N_CNTs, model, cnt_struct, cnt_coords, str_matrix, str_host):
	
	#Point at the center of the matrix
	Pc = (P0[0] + 0.5*Lxyz[0], P0[1] + 0.5*Lxyz[1], P0[2] + 0.5*Lxyz[2]);
	
	#Set for the matrix
	mdb.models[model].rootAssembly.Set(
		cells=mdb.models[model].rootAssembly.instances[str_matrix + '-1'].cells.findAt((Pc, ) ),
		name=str_host)
		
	#Variable to keep the count of CNT points
	acc_pts = 0
	
	#Sets for the CNTs
	for cnt_i in range(1, N_CNTs+1):
	
		#Get the string for the CNT part
		cnt_str = string_part('CNT', cnt_i)
		
		#Number of points in CNTi and its radius
		N_p = int(cnt_struct[cnt_i][0])
		
		#Generate name of set for cnt_i
		set_str = cnt_str + '_EESet'
		
		#Create the set for cnt_i
		mdb.models[model].rootAssembly.Set(
			cells=mdb.models[model].rootAssembly.instances[cnt_str + '-1'].cells.findAt( (cnt_coords[acc_pts+1],) ),
			name=set_str)
		
		#Increase the number of accumulated points
		acc_pts += N_p

#This functions creates the constraints for the embedded elements
def embedded_elements_constraints(N_CNTs, model, str_matrix, str_host):
	
	#Iterate over the CNTs
	for cnt_i in range(1, N_CNTs+1):
	
		#Get the string for the CNT part
		cnt_str = string_part('CNT', cnt_i)
		
		#Get the string for the CNT set
		set_str = cnt_str + '_EESet'
	
		#For cnt_i create an embedded element constraint
		#mdb.models['Model-1'].EmbeddedRegion(
		#	absoluteTolerance=0.0, fractionalTolerance=0.05, toleranceMethod=BOTH, weightFactorTolerance=1e-06,
		#	embeddedRegion=
		#		Region(cells=mdb.models['Model-1'].rootAssembly.instances[cnt_str + '-1'].cells.getSequenceFromMask(mask=('[#1 ]', ), )), 
		#	hostRegion=
		#		Region(cells=mdb.models['Model-1'].rootAssembly.instances[str_matrix + '-1'].cells.getSequenceFromMask(mask=('[#1 ]', ), )),
		#	name='EE-Constraint-%d' %(cnt_i))
		mdb.models[model].EmbeddedRegion(
			absoluteTolerance=0.0, fractionalTolerance=0.05, toleranceMethod=BOTH, weightFactorTolerance=1e-06,
			embeddedRegion=mdb.models[model].rootAssembly.sets[set_str],
			hostRegion=mdb.models[model].rootAssembly.sets[str_host],
			name='EE-Constraint-%d' %(cnt_i))
	
#This function generates the mesh for the matrix
def generate_matrix_mesh(model, str_matrix, matrixMeshSize):
	
	#Mesh the matrix
	#deviationFactor and minSizeFactor have the default values from Abaqus
	mdb.models[model].parts[str_matrix].seedPart(deviationFactor=0.1, minSizeFactor=0.1, size=matrixMeshSize)
	mdb.models[model].parts[str_matrix].generateMesh()

#This function generates the mesh for the CNTs
def generate_cnt_meshes(N_CNTs, cnt_struct):
	
	#Go through every CNT to mesh each of them
	for cnt_i in range(1, N_CNTs+1):
	
		#Get the string for the CNT part
		cnt_str = string_part('CNT', cnt_i)
		
		#Mesh cnt_i, use its radius as the element size
		#deviationFactor and minSizeFactor have the default values from Abaqus
		mdb.models['Model-1'].parts[cnt_str].seedPart(deviationFactor=0.1, minSizeFactor=0.1, size=cnt_struct[cnt_i][1])
		mdb.models['Model-1'].parts[cnt_str].generateMesh()
		
	#This seems to be required by Abaqus
	mdb.models['Model-1'].rootAssembly.regenerate()

#This function generates the mesh for the GSs
def generate_gs_meshes(N_GSs, model, code, meshSize):
    
    #Go through every CNT to mesh each of them
    for i in range(N_GSs):
    
        #Get the string for the CNT part
        str = string_part('GS',i)
        
        #Mesh GS i
        Mesh_GS(model, str, code, meshSize)
        
    #This seems to be required by Abaqus
    mdb.models[model].rootAssembly.regenerate()

#Mesh graphene sheets
def Mesh_GS(model, str_part, code, meshSize):
    mdb.models[model].parts[str_part].setElementType(
    elemTypes=(
        ElemType(elemCode=code[0], elemLibrary=STANDARD, secondOrderAccuracy=OFF, distortionControl=DEFAULT), 
        ElemType(elemCode=code[1], elemLibrary=STANDARD), 
        ElemType(elemCode=code[2], elemLibrary=STANDARD)), 
    regions=(mdb.models[model].parts[str_part].cells.getSequenceFromMask(('[#1 ]', ), ), ))

    mdb.models[model].parts[str_part].seedPart(deviationFactor=0.1, minSizeFactor=0.1, size=meshSize)
    mdb.models[model].parts[str_part].generateMesh()
    
#Create the embedded elements constrains for GSs
def create_embedded_elments_gs(N_GSs, model, matrixName, filler):
    
    #Iterate over the GSs
    for i in range(N_GSs):
        
        ##Create the constraint 
        Create_Embedded_Elements(model, matrixName, filler, i)
        
#Create embedded elements
def Create_Embedded_Elements(model, matrixName, filler, i):
        
    #Get the instance name
    str_inst = string_instance(filler, i)
        
    #Get the string for the embedded element constraint for GS i
    str_ee = ee_string(filler, i)
    
    #Add the embedded element constraint
    mdb.models[model].EmbeddedRegion(
        absoluteTolerance=0.0, 
        embeddedRegion=Region(
            cells=mdb.models[model].rootAssembly.instances[str_inst].cells.getSequenceFromMask(mask=('[#1 ]', ), )), 
        fractionalTolerance=0.05, 
        hostRegion=Region(
            cells=mdb.models[model].rootAssembly.instances[matrixName + '-1'].cells.getSequenceFromMask(mask=('[#1 ]', ), )), 
        name=str_ee, 
        toleranceMethod=BOTH, 
        weightFactorTolerance=1e-06)

def create_step_and_pbcs(model, P0, Lxyz, str_matrix, disp):
    
    #Create a step with default values
    #mdb.models['Model-1'].StaticStep(name='Step-1', previous='Initial', initialInc=0.1)
    #Create a step with fixed step
    mdb.models[model].StaticStep(
        initialInc=0.2, name='Step-1', noStop=OFF, 
        previous='Initial', timeIncrementationMethod=FIXED)
    
    #Calculate maximum coordinates of sample
    xmax = P0[0] + Lxyz[0]
    ymax = P0[1] + Lxyz[1]
    zmax = P0[2] + Lxyz[2]
    
    #Calculate midpoints
    xmid = P0[0] + 0.5*Lxyz[0]
    ymid = P0[1] + 0.5*Lxyz[1]
    zmid = P0[2] + 0.5*Lxyz[2]
    
    # Create reference points, which shuld be located outside the sample
    RPx = mdb.models[model].rootAssembly.ReferencePoint(point=(1.5*xmax, ymid, zmid))
    RPy = mdb.models[model].rootAssembly.ReferencePoint(point=(xmid, 1.5*ymax, zmid))
    RPz = mdb.models[model].rootAssembly.ReferencePoint(point=(xmid, ymid, 1.5*zmax))
    
    #Create sets for the reference points
    mdb.models[model].rootAssembly.Set(
        name='RP_X', 
        referencePoints=(mdb.models['Model-1'].rootAssembly.referencePoints[RPx.id], ))
    mdb.models[model].rootAssembly.Set(
        name='RP_Y', 
        referencePoints=(mdb.models['Model-1'].rootAssembly.referencePoints[RPy.id], ))
    mdb.models[model].rootAssembly.Set(
        name='RP_Z', 
        referencePoints=(mdb.models['Model-1'].rootAssembly.referencePoints[RPz.id], ))
    
    #Create sets for each face of the sample
    mdb.models[model].rootAssembly.Set(
        faces=mdb.models['Model-1'].rootAssembly.instances['Matrix-1'].faces.findAt(((xmax, ymid, zmid), )), 
        name='F_Xp')
    mdb.models[model].rootAssembly.Set(
        faces=mdb.models['Model-1'].rootAssembly.instances['Matrix-1'].faces.findAt(((P0[0], ymid, zmid), )), 
        name='F_Xn')
    mdb.models[model].rootAssembly.Set(
        faces=mdb.models['Model-1'].rootAssembly.instances['Matrix-1'].faces.findAt(((xmid, ymax, zmid), )), 
        name='F_Yp')
    mdb.models[model].rootAssembly.Set(
        faces=mdb.models['Model-1'].rootAssembly.instances['Matrix-1'].faces.findAt(((xmid, P0[1], zmid), )), 
        name='F_Yn')
    mdb.models[model].rootAssembly.Set(
        faces=mdb.models['Model-1'].rootAssembly.instances['Matrix-1'].faces.findAt(((xmid, ymid, zmax), )), 
        name='F_Zp')
    mdb.models[model].rootAssembly.Set(
        faces=mdb.models['Model-1'].rootAssembly.instances['Matrix-1'].faces.findAt(((xmid, ymid, P0[2]), )), 
        name='F_Zn')
    
    #Add the equations that will result in preiodic boundary conditions 
    mdb.models[model].Equation(
        name='Xp', terms=(( 1.0, 'F_Xp', 1), (0.5, 'RP_X', 1)))
    mdb.models[model].Equation(
        name='Xn', terms=((-1.0, 'F_Xn', 1), (0.5, 'RP_X', 1)))
    mdb.models[model].Equation(
        name='Yp', terms=(( 1.0, 'F_Yp', 2), (0.5, 'RP_Y', 2)))
    mdb.models[model].Equation(
        name='Yn', terms=((-1.0, 'F_Yn', 2), (0.5, 'RP_Y', 2)))
    mdb.models[model].Equation(
        name='Zp', terms=(( 1.0, 'F_Zp', 3), (0.5, 'RP_Z', 3)))
    mdb.models[model].Equation(
        name='Zn', terms=((-1.0, 'F_Zn', 3), (0.5, 'RP_Z', 3)))
    
    #Create a table to obtain different steps of the simulation
    #mdb.models['Model-1'].TabularAmplitude(
    #    data=((0.0, 0.0), (0.1, 0.1), (0.2, 0.2), (0.3, 0.3), (0.4, 0.4), (0.5, 0.5), 
    #        (0.6, 0.6), (0.7, 0.7), (0.8, 0.8) , (0.9, 0.9), (1.0, 1.0)), 
    #    name='Amp-1', 
    #    smooth=SOLVER_DEFAULT, timeSpan=STEP)
    #Create a simple amplitude
    mdb.models[model].TabularAmplitude(
        data=((0.0, 0.0), (1.0, 1.0)), 
        name='Amp-1', 
        smooth=SOLVER_DEFAULT, timeSpan=STEP)
    
    #Apply the displacement BC to 
    mdb.models[model].DisplacementBC(
        amplitude='Amp-1', 
        createStepName='Step-1', 
        distributionType=UNIFORM, 
        fieldName='', fixed=OFF, localCsys=None, 
        name='BC-2', 
        region=Region(referencePoints=(mdb.models['Model-1'].rootAssembly.referencePoints[RPx.id], )), 
        u1=disp, u2=UNSET, u3=UNSET, ur1=UNSET, ur2=UNSET, ur3=UNSET)
#This function partitions the CNT into various cells
def datum_points(cnt_rad, cnt_start, cnt_end, cnt_coords, str_part):
	
	#Calculate the radius multiplied by cos(PI/4)
	new_rad = cnt_rad*cos45
	
	#Iterate over the CNT points, starting on the secont point and finishing on the previous to last point
	for i in range(cnt_start+1, cnt_end): 
		
		#print('i=',i-cnt_start,' cells=', len(mdb.models['Model-1'].parts[str_part].cells))
		
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
		
		#Get a vector that is the average of v1 and v2
		vm = ((v1[0]+v2[0])/2, (v1[1]+v2[1])/2, (v1[2]+v2[2])/2)
		#print('vm')
		
		#Get the normal for the plane at the midpoint
		N = cross(cross(v1, vm), vm)
		
		#Get the rotation matrix corresponding to the current CNT segment (P1P2)
		R = rotation_matrix((0.0, 0.0, 0.0), N)
		#print('R')
		
		#Calculate points on the octagon
		#(0.0, cnt_rad, 0.0)
		Q1 = (cnt_rad*R[0][1] + P2[0], cnt_rad*R[1][1] + P2[1], cnt_rad*R[2][1] + P2[2])
		#(new_rad, new_rad)
		Q2 = (new_rad*(R[0][0] + R[0][1]) + P2[0], new_rad*(R[1][0] + R[1][1]) + P2[1], new_rad*(R[2][0] + R[2][1]) + P2[2])
		#(cnt_rad, 0.0
		Q3 = (cnt_rad*R[0][0] + P2[0], cnt_rad*R[1][0] + P2[1], cnt_rad*R[2][0] + P2[2])
		#Create datum points
		mdb.models['Model-1'].parts[str_part].DatumPointByCoordinate(coords=Q1)
		mdb.models['Model-1'].parts[str_part].DatumPointByCoordinate(coords=Q2)
		mdb.models['Model-1'].parts[str_part].DatumPointByCoordinate(coords=Q3)
		mdb.models['Model-1'].parts[str_part].DatumPointByCoordinate(coords=P2)
		#print('P2=',P2)

		
#This function creates two sets for two of the vertices of the matrix (sample), where each set has only one node
#These nodes are on the diagonal of the cuboid that defines the matrix (sample) 
def create_sets_for_matrix(P0, Lxyz, str_matrix):
	
	#Create the set of the lower left corner
	#NOTE: The parenthesis, comma and the space in the nodes argument is beacuse a tuple is needed but
	#the operation getClosest returns an object. The parenthesis, comma and space make it a tuple
	mdb.models['Model-1'].rootAssembly.Set(
		nodes=mdb.models['Model-1'].rootAssembly.instances[str_matrix+"-1"].nodes.getByBoundingSphere(center=P0, radius=0.001),
		name='Matrix0')
	#This way does not allow to get displacements for some reason
	#I guess the mesh is on instance and not in part
	#mdb.models['Model-1'].parts[str_matrix].Set(
	#	nodes=(mdb.models['Model-1'].parts[str_matrix].nodes.getByBoundingSphere(center=P0, radius=0.001), ),
	#	name='Matrix0')
	
	#print('P0=',P0,' str_matrix=',str_matrix)
	#nd = mdb.models['Model-1'].rootAssembly.instances[str_matrix+"-1"].nodes.getByBoundingSphere(center=P0, radius=0.001)
	#print('nd0=',nd[0].coordinates)
	
	#Create opposite corner
	corner = (P0[0]+Lxyz[0], P0[1]+Lxyz[1], P0[2]+Lxyz[2])
	
	#Create the set of the opposite corner
	mdb.models['Model-1'].rootAssembly.Set(
		nodes=mdb.models['Model-1'].rootAssembly.instances[str_matrix+"-1"].nodes.getByBoundingSphere(center=corner, radius=0.001),
		name='Matrix1')
	#nd = mdb.models['Model-1'].rootAssembly.instances[str_matrix+"-1"].nodes.getByBoundingSphere(center=corner, radius=0.001)
	#print('nd1=',nd[0].coordinates)
		
	#Perform union of two sets
	#mdb.models['Model-1'].rootAssembly.SetByBoolean(
	#	name='Matrix2', 
	#	operation=UNION, 
	#	sets=(mdb.models['Model-1'].rootAssembly.sets['Matrix0'], mdb.models['Model-1'].rootAssembly.sets['Matrix1'])
	#	)

#This function creates a set for the nodes that correspond to the centerline of that CNT
def create_set_for_cnt_points(cnt_i, cnt_rad, cnt_start, cnt_end, cnt_coords, cnt_str):
	
	#Calculate radius for serching nodes
	new_rad = cnt_rad/4;
	
	#Get the name of the node set
	node_set_str = string_node_set('CNT', cnt_i)
	
	#Create a set with the name of the node set and containing the first point in the CNT
	mdb.models['Model-1'].rootAssembly.Set(
		nodes=mdb.models['Model-1'].rootAssembly.instances[cnt_str].nodes.getByBoundingSphere(center=cnt_coords[cnt_start], radius=new_rad),
		name=node_set_str)
	
	#Iterate over all points of the CNT, starting on the second point since the first point is already in the set
	for i in range(cnt_start+1, cnt_end+1):
		
		#Create a temporary set with node i
		node_set_tmp = mdb.models['Model-1'].rootAssembly.Set(
			nodes=mdb.models['Model-1'].rootAssembly.instances[cnt_str].nodes.getByBoundingSphere(center=cnt_coords[i], radius=new_rad),
			name='tmp_set')
		
		#Merge two sets
		#mdb.models['Model-1'].rootAssembly.SetByMerge(
		#	name=node_set_str, 
		#	sets=(node_set_tmp, mdb.models['Model-1'].rootAssembly.sets[node_set_str]))
		
		#Perform union of two sets
		mdb.models['Model-1'].rootAssembly.SetByBoolean(
			name=node_set_str, 
			operation=UNION, 
			sets=(node_set_tmp, mdb.models['Model-1'].rootAssembly.sets[node_set_str]))
		
	#Print the length of the set
	print('%s nodes=%d points=%d'%(node_set_str, len(mdb.models['Model-1'].rootAssembly.sets[node_set_str].nodes), cnt_end+1-cnt_start))

#This function creates all sets for the nodes that correspond to the centerline of that CNT
def create_all_sets_for_cnt_points(N_CNTs, cnt_struct, cnt_coords):
	
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
		create_set_for_cnt_points(cnt_i, cnt_rad, acc_pts, acc_pts+N_p-1, cnt_coords, cnt_str)
		
		#Increase the number of accumulated points
		acc_pts += N_p
	
	
######################################---ABAQUS FUNCTIONS---########################################
######################################----------END---------########################################



######################################----------END---------########################################



######################################---MAIN PROGRAM---########################################


#########################################---OPENING CSV FILES---#########################################
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
    
#Get the number of CNTs
N_CNTs = int(cnt_struct[0][0])
print('The number of CNTs is:', N_CNTs)


with open(gs_file) as f:
    #x = csv.reader(f, quoting=csv.QUOTE_NONNUMERIC)
    gs_data = list(csv.reader(f, quoting=csv.QUOTE_NONNUMERIC))

#Get the number of GSs
N_GSs = len(gs_data)
print('The number of GSs is:', N_GSs)

#Read the sample geometry
with open(sample_file) as f:
    (P0, Lxyz) = list(csv.reader(f, quoting=csv.QUOTE_NONNUMERIC))
#print(P0)
#print(Lxyz)


#######################################---MODEL CREATION---#########################################


###----------------------------------Start: Create Parts----------------------------------###

#Generate the matrix part
matrix_part(P0, Lxyz, str_matrix, sheetSize)

#A small fixed number of CNTs used for debugging and testing
#N_CNTs=2

if N_CNTs != 0:
	#Generate all CNT parts
	cnt_parts_all(N_CNTs, cnt_struct, cnt_coords)

#Box used for cutting GSs
#Create_BiggerBox(modelName, sheetSize, biggerBoxName, matrixLength)

if N_GSs != 0:
	#Create grpahene sheet parts
	gs_parts_all(N_GSs, gs_data, modelName, sheetSize)
    
###----------------------------------End: Create Parts----------------------------------###


###----------------------------------Start: Materials and sections----------------------------------###
#Create matrix material
#Create_Material(modelName, matrixMaterial, matrixDensity, matrixModulus, matrixPoissonR, matrixExpCoeff, matrixThermConductivity, matrixSpecHeat, matrixElecConductivity)
#Create GS material
#Create_Material(modelName, fillerMaterial, fillerDensity, fillerModulus, fillerPoissonR, fillerExpCoeff, fillerThermConductivity, fillerSpecHeat, fillerElecConductivity)

#Creating a section for matrix
#Create_Section(modelName, matrixMaterial)
#Create a section for GS
#Create_Section(modelName, fillerMaterial)

#Generate materials and assign sections
materials_and_sections_matrix(modelName, str_matrix, matrixMaterial, matrixSection, matrixDensity, matrixModulus, matrixPoissonR)

#Create material and section for CNTs (only elastic properties)
materials_and_sections_cnt(modelName, N_CNTs, cntMaterial, cntSection, cntDensity, cntModulus, cntPoissonR)

#Create material and section for GS (only elastic properties)
materials_and_sections_gs(modelName, N_GSs, gsMaterial, gsSection, gsDensity, gsModulus, gsPoissonR)


###----------------------------------End: Materials and sections----------------------------------###


###----------------------------------Start: Create assembly----------------------------------###

#Generate assembly
generate_assembly(N_CNTs, N_GSs, modelName, str_matrix, gs_data)

###-----------------------------------End: Create assembly-----------------------------------###

###-------------------------------------Start: BCs-------------------------------------###

#Defining boundary conditions
#mdb.models[modelName].CoupledThermalElectricalStructuralStep(deltmx=0.1, initialInc=0.01, name='Step_name', previous='Initial')

#Boundary conditions for testing
# Displacement_BoundaryConditions(modelName, matrixName)
# Temperature_BoundaryConditions(modelName, matrixName, tempApplied, numRows)

# for numPart in indexInside:
#     #Create embedded elements (GS) one by one into the host region (matrix)
#     Create_Embedded_Elements(modelName, matrixName, numPart)
#     #Mesh each graphene sheet (meshing: sizeGS < sizeMatrix)
#     Mesh_GS(modelName, selectedElementCode, eeMeshSize, numPart)
#
# for numPart in indexOutside:
#     #Create embedded elements (GS) one by one into the host region (matrix)
#     Create_Embedded_Elements_CuttedGS(modelName, matrixName, numPart)
#     #Mesh each graphene sheet (meshing: sizeGS < sizeMatrix)
#     Mesh_GS(modelName, selectedElementCode, eeMeshSize, numPart)

if N_GSs != 0:
	#Create all embedded elements constraints for GSs when using option all_in
	create_embedded_elments_gs(N_GSs, modelName, str_matrix, 'GS')

if N_CNTs != 0:
	#Create sets that will be used when creating the embedded element constraints for CNTs
	sets_for_embedded_elements(P0, Lxyz, N_CNTs, modelName, cnt_struct, cnt_coords, str_matrix, str_host)
	
	#Create embedded elements constraints for CNTs
	embedded_elements_constraints(N_CNTs, modelName, str_matrix, str_host)

#Create Step and add boundary conditions
create_step_and_pbcs(modelName, P0, Lxyz, str_matrix, disp)

###-------------------------------------End: BCs-------------------------------------###

###-------------------------------------Start: Meshing-------------------------------------###

#Calculate the matrix mesh size

#Generate meshe for the matrix
generate_matrix_mesh(modelName, str_matrix, matrixMeshSize)

if N_CNTs != 0:
	
	#Generate mesh for CNTs
	generate_cnt_meshes(N_CNTs, cnt_struct)
	
	#Create sets for central CNT nodes
	#NOTE: Sets are generated on root assembly
	create_all_sets_for_cnt_points(N_CNTs, cnt_struct, cnt_coords)

if N_GSs != 0:
	
	#Generate mesh for GSs when using option all_in
	generate_gs_meshes(N_GSs, modelName, selectedElementCode, eeMeshSize)

#Create sets for the matrix (sample in the C++ code)
#NOTE: Sets are generated on root assembly
create_sets_for_matrix(P0, Lxyz, str_matrix)



#Mesh the matrix (its size should be bigger than the meshing size of the GS). Meshes are being created with code
# mdb.models[modelName].parts[matrixName].setElementType(
#     elemTypes=(
#         ElemType(elemCode=selectedElementCode[0], elemLibrary=STANDARD, secondOrderAccuracy=OFF, distortionControl=DEFAULT), 
#         ElemType(elemCode=selectedElementCode[1], elemLibrary=STANDARD), 
#         ElemType(elemCode=selectedElementCode[2], elemLibrary=STANDARD)), 
#     regions=(mdb.models[modelName].parts[matrixName].cells.getSequenceFromMask(('[#1 ]', ), ), ))
# mdb.models[modelName].parts[matrixName].seedPart(deviationFactor=0.1, minSizeFactor=0.1, size=matrixMeshSize)
# mdb.models[modelName].parts[matrixName].generateMesh()
# mdb.models[modelName].rootAssembly.regenerate()
###-------------------------------------End: Meshing--------------------------------------###


#Create and submit job using Abaqus default values
mdb.Job(atTime=None, contactPrint=OFF, description='', echoPrint=OFF, 
	explicitPrecision=SINGLE, getMemoryFromAnalysis=True, historyPrint=OFF, 
	memory=90, memoryUnits=PERCENTAGE, model='Model-1', modelPrint=OFF, 
	multiprocessingMode=DEFAULT, name='Job-1', nodalOutputPrecision=SINGLE, 
	numCpus=1, numGPUs=0, queue=None, resultsFormat=ODB, scratch='', type=
	ANALYSIS, userSubroutine='', waitHours=0, waitMinutes=0)
mdb.jobs['Job-1'].submit(consistencyChecking=OFF)