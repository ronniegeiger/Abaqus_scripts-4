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

######################################---ABAQUS FUNCTIONS---########################################

def matrix_part(P0, Lxyz, str_matrix):
	
	#Create a sketch
	mdb.models['Model-1'].ConstrainedSketch(name='__profile__', sheetSize=200.0)
	
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
    		    
#This function creates a CNT in Part module
def cnt_part(cnt_i, cnt_rad, cnt_start, cnt_end, cnt_coords):
	
	#Get the string for the CNT
	str_part = cnt_string_part(cnt_i)
	
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

#This function creates the assembly by creating the instances of each part
def generate_assembly(N_CNTs, str_matrix):
	
	#Necessary command to create assembly
	mdb.models['Model-1'].rootAssembly.DatumCsysByDefault(CARTESIAN)
	
	#Generate name for matrix instance
	str_mat_inst = '%s-1' % (str_matrix)
	
	#Create instance for matrix
	mdb.models['Model-1'].rootAssembly.Instance(dependent=ON, name=str_mat_inst,
		part=mdb.models['Model-1'].parts[str_matrix])
	
	#Create instances for each CNT
	for i in range(1, N_CNTs+1):
		
		#Generate name for CNT part
		str_cnt = cnt_string_part(i)
		
		#Generate name for CNT instance
		str_cnt_inst = cnt_string_instance(i)
		
		#Generate CNT instance
		mdb.models['Model-1'].rootAssembly.Instance(dependent=ON, name=str_cnt_inst,
			part=mdb.models['Model-1'].parts[str_cnt])		

#This function creates the matrix material and assigns it to a section
def materials_and_sections_matrix(str_matrix, matrixMaterial, matrixSection, matrixDensity, matrixModulus, matrixPoissonR):
	
	#Create matrix material
	mdb.models['Model-1'].Material(description='Polymer', name=matrixMaterial)
	mdb.models['Model-1'].materials[matrixMaterial].Elastic(table=((matrixModulus, matrixPoissonR), ))
	
	#Assign material to section
	mdb.models['Model-1'].HomogeneousSolidSection(material=matrixMaterial, name=matrixSection, thickness=None)
	
	#Assign section to matrix
	#print('cells=',len(mdb.models['Model-1'].parts[str_matrix].cells))
	mdb.models['Model-1'].parts[str_matrix].SectionAssignment(offset=0.0, offsetField='', offsetType=MIDDLE_SURFACE, 
		region=Region(cells=mdb.models['Model-1'].parts[str_matrix].cells), 
		sectionName=matrixSection, thicknessAssignment=FROM_SECTION)

#This function creates the CNT material and assigns it to a section
def materials_and_sections_cnt(N_CNTs, cntMaterial, cntSection, cntDensity, cntModulus, cntPoissonR):
	
	#Create CNT material
	mdb.models['Model-1'].Material(name=cntMaterial)
	mdb.models['Model-1'].materials[cntMaterial].Elastic(table=((cntModulus, cntPoissonR), ))
	
	#Assign material to section
	mdb.models['Model-1'].HomogeneousSolidSection(material=cntMaterial, name=cntSection, thickness=None)
	
	#Iterate over the CNTs
	for cnt_i in range(1, N_CNTs+1):
	
		#Get the string for the CNT part
		cnt_str = cnt_string_part(cnt_i)
		
		#Assign the CNT section to cnt_i
		mdb.models['Model-1'].parts[cnt_str].SectionAssignment(offset=0.0, offsetField='', offsetType=MIDDLE_SURFACE, 
			region=Region(cells=mdb.models['Model-1'].parts[cnt_str].cells.getSequenceFromMask(mask=('[#1 ]', ), )),
			sectionName=cntSection, thicknessAssignment=FROM_SECTION)

#This function creates the sets that will be used when creating the embedded element constraints
def sets_for_embedded_elements(P0, Lxyz, N_CNTs, cnt_struct, cnt_coords, str_matrix, str_host):
	
	#Point at the center of the matrix
	Pc = (P0[0] + 0.5*Lxyz[0], P0[1] + 0.5*Lxyz[1], P0[2] + 0.5*Lxyz[2]);
	
	#Set for the matrix
	mdb.models['Model-1'].rootAssembly.Set(
		cells=mdb.models['Model-1'].rootAssembly.instances[str_matrix + '-1'].cells.findAt((Pc, ) ),
		name=str_host)
		
	#Variable to keep the count of CNT points
	acc_pts = 0
	
	#Sets for the CNTs
	for cnt_i in range(1, N_CNTs+1):
	
		#Get the string for the CNT part
		cnt_str = cnt_string_part(cnt_i)
		
		#Number of points in CNTi and its radius
		N_p = int(cnt_struct[cnt_i][0])
		
		#Generate name of set for cnt_i
		set_str = cnt_str + '_EESet'
		
		#Create the set for cnt_i
		mdb.models['Model-1'].rootAssembly.Set(
			cells=mdb.models['Model-1'].rootAssembly.instances[cnt_str + '-1'].cells.findAt( (cnt_coords[acc_pts+1],) ),
			name=set_str)
		
		#Increase the number of accumulated points
		acc_pts += N_p

#This functions creates the constraints for the embedded elements
def embedded_elements_constraints(N_CNTs, str_matrix, str_host):
	
	#Iterate over the CNTs
	for cnt_i in range(1, N_CNTs+1):
	
		#Get the string for the CNT part
		cnt_str = cnt_string_part(cnt_i)
	
		#Get the string for the CNT part
		cnt_str = cnt_string_part(cnt_i)
		
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
		mdb.models['Model-1'].EmbeddedRegion(
			absoluteTolerance=0.0, fractionalTolerance=0.05, toleranceMethod=BOTH, weightFactorTolerance=1e-06,
			embeddedRegion=mdb.models['Model-1'].rootAssembly.sets[set_str],
			hostRegion=mdb.models['Model-1'].rootAssembly.sets[str_host],
			name='EE-Constraint-%d' %(cnt_i))
	
#This function generates the mesh for the matrix
def generate_matrix_mesh(str_matrix, Lxyz, cnt_rad_max):
	
	#Determine the size of the mesh for the polymer matrix
	#First determine the size of the elements considering 10 elements per side
	#In case the sample is not a cube, then first I need to find the shortest side
	#This side will determine the size of the element
	el_size_mat = min(Lxyz)/10.0
	
	#Check that the element size is at least twice the maximum CNT radius
	if (el_size_mat - 2*cnt_rad_max <= Zero ):
		
		#The element size for the matrix is too small, so set it to be twice the maximum CNT radius
		el_size_mat = 2*cnt_rad_max
	
	#Mesh the matrix
	#deviationFactor and minSizeFactor have the default values from Abaqus
	mdb.models['Model-1'].parts[str_matrix].seedPart(deviationFactor=0.1, minSizeFactor=0.1, size=el_size_mat)
	mdb.models['Model-1'].parts[str_matrix].generateMesh()

#This function generates the mesh for the CNTs
def generate_cnt_meshes(N_CNTs, cnt_struct):
	
	#Go through every CNT to mesh each of them
	for cnt_i in range(1, N_CNTs+1):
	
		#Get the string for the CNT part
		cnt_str = cnt_string_part(cnt_i)
		
		#Mesh cnt_i, use its radius as the element size
		#deviationFactor and minSizeFactor have the default values from Abaqus
		mdb.models['Model-1'].parts[cnt_str].seedPart(deviationFactor=0.1, minSizeFactor=0.1, size=cnt_struct[cnt_i][1])
		mdb.models['Model-1'].parts[cnt_str].generateMesh()
		
	#This seems to be required by Abaqus
	mdb.models['Model-1'].rootAssembly.regenerate()
	
#This functions creates a step and adds the boundary conditions for tension along the z-axis
def create_step_and_bcs(str_matrix, disp):
	
	#Create a step with default values
	mdb.models['Model-1'].StaticStep(name='Step-1', previous='Initial')
	
	#Add the displacement boundary condition to the "top" face of the sample
	mdb.models['Model-1'].rootAssembly.Set(
		faces=mdb.models['Model-1'].rootAssembly.instances[str_matrix + '-1'].faces.getSequenceFromMask( ('[#10 ]', ), ), 
		name='Dips_BC')
	mdb.models['Model-1'].DisplacementBC(
		amplitude=UNSET, createStepName='Step-1', distributionType=UNIFORM, fieldName='', fixed=OFF, localCsys=None, 
		name='BC-1', 
		region=mdb.models['Model-1'].rootAssembly.sets['Dips_BC'], 
		u1=UNSET,u2=UNSET, u3=disp, ur1=UNSET, ur2=UNSET, ur3=UNSET)
	
	#Add the pinned boundary condition to the "bottom" face of the sample
	mdb.models['Model-1'].rootAssembly.Set(
		faces=mdb.models['Model-1'].rootAssembly.instances[str_matrix + '-1'].faces.getSequenceFromMask(('[#20 ]', ), ),
		name='Encastre_BC')
	mdb.models['Model-1'].EncastreBC(createStepName='Step-1', localCsys=None, name='BC-2', 
		region=mdb.models['Model-1'].rootAssembly.sets['Encastre_BC'])


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
	node_set_str = cnt_string_node_set(cnt_i)
	
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
		cnt_str = cnt_string_instance(cnt_i)
		
		#Create a set for the nodes at the centerline of cnt_i
		create_set_for_cnt_points(cnt_i, cnt_rad, acc_pts, acc_pts+N_p-1, cnt_coords, cnt_str)
		
		#Increase the number of accumulated points
		acc_pts += N_p
	
	
######################################---ABAQUS FUNCTIONS---########################################
######################################----------END---------########################################



######################################---GLOBAL VARIABLES---########################################

#Files to read
coord_file = 'C:\Users\grafenicas\Documents\CNTs_embedde_elements\Networks\cnt_coordinates.csv'
struct_file = 'C:\Users\grafenicas\Documents\CNTs_embedde_elements\Networks\cnt_struct.csv'
sample_file = 'C:\Users\grafenicas\Documents\CNTs_embedde_elements\Networks\sample_geom.csv'
	
#Define the string for the name of the matrix part
str_matrix = 'Matrix'

#String for the host set (i.e., the matrix)
str_host = 'host_Set'

#Cosine of 45 degrees (PI/4)
cos45 = 0.7071067812

#"Zero" for comparing floating point numbers
Zero = 1e-7

#Displacement to be applied (microns)
disp = 0.5

#Matrix properties
#Define the name and section of the matrix
matrixMaterial = 'Matrix_mat'
matrixSection = 'Matrix_sec'
#Define the mass density (kg/m3)
matrixDensity = 905
#Define the elastic modulus (GPa)
matrixModulus = 1e6
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
cntMaterial = 'CNT_mat'
cntSection = 'CNT_sec'
#Define the mass density (kg/m3)
cntDensity = 2200
#Define the elastic modulus (GPa)
cntModulus = 1e9
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
cnt_rad_max = 0.005

######################################---GLOBAL VARIABLES---########################################
######################################----------END---------########################################



######################################---MAIN PROGRAM---########################################

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

#Generate the matrix part
matrix_part(P0, Lxyz, str_matrix)
    
#Get the number of CNTs
N_CNTs = int(cnt_struct[0][0])

#A small fixed number used for debugging and testing
#N_CNTs=2

#Generate all CNT parts
cnt_parts_all(N_CNTs, cnt_struct, cnt_coords)

#Generate materials and assign sections
materials_and_sections_matrix(str_matrix, matrixMaterial, matrixSection, matrixDensity, matrixModulus, matrixPoissonR)
materials_and_sections_cnt(N_CNTs, cntMaterial, cntSection, cntDensity, cntModulus, cntPoissonR)

#Generate assembly
generate_assembly(N_CNTs, str_matrix)

#Create sets that will be used when creating the embedded element constraints
sets_for_embedded_elements(P0, Lxyz, N_CNTs, cnt_struct, cnt_coords, str_matrix, str_host)

#Create embedded elements constraints
embedded_elements_constraints(N_CNTs, str_matrix, str_host)

#Create Step and add boundary conditions
create_step_and_bcs(str_matrix, disp)

#Generate meshes
generate_matrix_mesh(str_matrix, Lxyz, cnt_rad_max)
generate_cnt_meshes(N_CNTs, cnt_struct)

#Create sets for the matrix (sample in the C++ code)
#NOTE: Sets are generated on root assembly
create_sets_for_matrix(P0, Lxyz, str_matrix)

#Create sets for central CNT nodes
#NOTE: Sets are generated on root assembly
create_all_sets_for_cnt_points(N_CNTs, cnt_struct, cnt_coords)

#Create and submit job using Abaqus default values
mdb.Job(atTime=None, contactPrint=OFF, description='', echoPrint=OFF, 
	explicitPrecision=SINGLE, getMemoryFromAnalysis=True, historyPrint=OFF, 
	memory=90, memoryUnits=PERCENTAGE, model='Model-1', modelPrint=OFF, 
	multiprocessingMode=DEFAULT, name='Job-1', nodalOutputPrecision=SINGLE, 
	numCpus=1, numGPUs=0, queue=None, resultsFormat=ODB, scratch='', type=
	ANALYSIS, userSubroutine='', waitHours=0, waitMinutes=0)
mdb.jobs['Job-1'].submit(consistencyChecking=OFF)