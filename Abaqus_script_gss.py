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

import math
import csv

##########################################---VARIABLES---##########################################

#Add the path of the .csv file
gs_file = 'C:\Users\grafenicas\Documents\CNTs_embedde_elements\Job-testing-GS\gnp_data.csv'
sample_file = 'C:\Users\grafenicas\Documents\CNTs_embedde_elements\Job-testing-GS\sample_geom.csv'
    
#Define the string for the name of the matrix part
str_matrix = 'Matrix'

#String for the host set (i.e., the matrix)
str_host = 'host_Set'

#Define the sheetsize for the drawing space
sheetSize = 20.0
#Mesh size for the matrix (um)
matrixMeshSize = 2
#Define the EmbeddedMesh/HostMesh ratio
meshRatio = 0.75
#Temperature difference (C)
tempApplied = 50.0

#Define the name of the matrix
matrixMaterial = 'Polypropylene'
matrixSection = 'Matrix_sec'
#Define the mass density (kg/m3)
matrixDensity = 905
#Define the elastic modulus (GPa)
matrixModulus = 1.3
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

#Define the name of the filler
gsMaterial = 'Graphene'
gsSection = 'Graphene-sec'
#Define the mass density (kg/m3)
gsDensity = 2200
#Define the elastic modulus (GPa)
gsModulus = 1000
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

#Displacement to be applied (microns)
disp = -0.15

#Analsys type for the element
# 1 = coupled temperature-displacement
# 2 = coupled thermal-electrical-structural
# 3 = thermal-electric
# 4 = only displacement
elemTypeCode = 4

##########################################---CONSTANT---##########################################
rad_to_deg = 180/math.pi

######################################---ABAQUS FUNCTIONS---########################################

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

#Create matrix
def matrix_part(P0, Lxyz, str_matrix, sheetsz):
    
    #Create a sketch
    mdb.models['Model-1'].ConstrainedSketch(name='__profile__', sheetSize=sheetsz)
    
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

def Create_BiggerBox(model, sheetsz, part, length):

    GS_geometry = data[1]
    
    #Get the length along each direction
    length_x = GS_geometry[0]
    length_y = GS_geometry[1]
    thickness = GS_geometry[2]

    margin = math.ceil(math.sqrt(length_x*length_x+length_y*length_y+thickness*thickness))
    lengthBiggerBox = length+2*margin

    #Create the matrix geometry
    mdb.models[model].ConstrainedSketch(name='__profile__', sheetSize=sheetsz)
    mdb.models[model].sketches['__profile__'].rectangle(point1=(0.0, 0.0), point2=(lengthBiggerBox, lengthBiggerBox))
    mdb.models[model].Part(dimensionality=THREE_D, name=part, type=DEFORMABLE_BODY)

    mdb.models[model].parts[part].BaseSolidExtrude(depth=lengthBiggerBox, sketch=mdb.models[model].sketches['__profile__'])
    del mdb.models[model].sketches['__profile__']

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

#This function creates the matrix material and assigns it to a section
def materials_and_sections_matrix(modelName, str_matrix, matrixMaterial, matrixSection, matrixDensity, matrixModulus, matrixPoissonR):
    
    #Create matrix material
    mdb.models[modelName].Material(description='Polymer', name=matrixMaterial)
    mdb.models[modelName].materials[matrixMaterial].Elastic(table=((matrixModulus, matrixPoissonR), ))
    
    #Assign material to section
    create_section(modelName, matrixMaterial, matrixSection)
    
    #Assign section to matrix
    #print('cells=',len(mdb.models['Model-1'].parts[str_matrix].cells))
    mdb.models[modelName].parts[str_matrix].SectionAssignment(offset=0.0, offsetField='', offsetType=MIDDLE_SURFACE, 
        region=Region(cells=mdb.models[modelName].parts[str_matrix].cells), 
        sectionName=matrixSection, thicknessAssignment=FROM_SECTION)

#This function creates the GS material and assigns it to a section
def materials_and_sections_gs(N_GSs, gsMaterial, gsSection, gsDensity, gsModulus, gsPoissonR):
    
    #Create GS material
    mdb.models['Model-1'].Material(name=gsMaterial)
    mdb.models['Model-1'].materials[gsMaterial].Elastic(table=((gsModulus, gsPoissonR), ))
    
    #Assign material to section
    create_section('Model-1', gsMaterial, gsSection)
    
    #Iterate over the CNTs
    for gs_i in range(N_GSs):
    
        #Get the string for the GS part
        gs_str = string_part('GS', gs_i)
        
        #Assign the CNT section to gs_i
        mdb.models['Model-1'].parts[gs_str].SectionAssignment(
            offset=0.0, offsetField='', offsetType=MIDDLE_SURFACE, 
            region=Region(cells=mdb.models['Model-1'].parts[gs_str].cells.getSequenceFromMask(mask=('[#1 ]', ), )),
            sectionName=gsSection, thicknessAssignment=FROM_SECTION)
        
#Create material and assign the properties
def Create_Material(model, materialName, massDensity, elasModulus, poissonRatio, expanCoefficient, thermConductivity, specHeat, elecConductivity):
    mdb.models[model].Material(name=materialName)
    mdb.models[model].materials[materialName].Density(table=((massDensity, ), ))
    mdb.models[model].materials[materialName].Elastic(table=((elasModulus*1e9, poissonRatio), ))
    mdb.models[model].materials[materialName].Expansion(table=((expanCoefficient*1e-5, ), ))
    mdb.models[model].materials[materialName].Conductivity(table=((thermConductivity, ), ))
    mdb.models[model].materials[materialName].SpecificHeat(table=((specHeat, ), ))
    mdb.models[model].materials[materialName].ElectricalConductivity(table=((elecConductivity, ), ))

#Create section
def create_section(model, materialName, sectionName):
    mdb.models[model].HomogeneousSolidSection(
        material=materialName, 
        name=sectionName, 
        thickness=None)
   
#This function generates all instances and the assembly  
def generate_assembly(N_GSs, model, str_matrix, gs_data):
    
    #Necessary command to create assembly
    mdb.models[model].rootAssembly.DatumCsysByDefault(CARTESIAN)
    
    #Generate name for matrix instance
    str_mat_inst = '%s-1' % (str_matrix)
    
    #Create instance for matrix
    mdb.models[model].rootAssembly.Instance(dependent=ON, name=str_mat_inst,
        part=mdb.models[model].parts[str_matrix]) 
    
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
    #     '[#20301 #0 #201000 #0 #3010000 #2 ]', ), ))

def Create_CuttingBox(model, partMatrix, partBox, GS_geometry):
    
    #Get the paramters for cutting a GSs (for readability of the script)
    length_x = GS_geometry[0]
    length_y = GS_geometry[1]
    thickness = GS_geometry[2]

    margin = math.ceil(math.sqrt(length_x*length_x+length_y*length_y+thickness*thickness))

    mdb.models[model].rootAssembly.Instance(dependent=OFF, name=partBox + '-1', 
        part=mdb.models[model].parts[partBox])
    mdb.models[model].rootAssembly.translate(instanceList=(partBox + '-1', ), 
        vector=(-margin, -margin, -margin))

    mdb.models[model].rootAssembly.InstanceFromBooleanCut(cuttingInstances=(
        mdb.models[model].rootAssembly.instances[partMatrix + '-1'], ), 
        instanceToBeCut=
        mdb.models[model].rootAssembly.instances[partBox + '-1'], name='CUTTER'
        , originalInstances=SUPPRESS)

    mdb.models[model].rootAssembly.features[partMatrix + '-1'].resume()

indexOutside = []
indexInside = []

def Cut_GS(model, i, str_inst, str_part):
    
    #Get the root assembly
    GS_assembly = mdb.models[model].rootAssembly
    
    #Get all vertices from a graphene sheet
    GS_vertex = GS_assembly.instances[str_inst].vertices
    
    #Initialize with false
    checkOutside = False
    
    for numVertex in range(0,8):
        #Get coordinates of a specific vertex from a graphene sheet
        GS_vertex_coord = GS_assembly.getCoordinates(GS_vertex[numVertex])

        for indexCoord in range(0,3):

            if GS_vertex_coord[indexCoord] > matrixLength or  GS_vertex_coord[indexCoord] < 0:
                checkOutside = True
                break
        
        if checkOutside == True:
            GS_assembly.InstanceFromBooleanCut(
                cuttingInstances=(GS_assembly.instances['CUTTER-1'], ), 
                instanceToBeCut=GS_assembly.instances[str_inst], name=str_part, 
                originalInstances=SUPPRESS)
            GS_assembly.features['CUTTER-1'].resume()

            indexOutside.append(i)

            #Go to next GS if one coordinate of a vertex is outside RVE
            break

    if checkOutside == False:
        indexInside.append(i)    

    if i == numRows - 1:
        GS_assembly.features['CUTTER-1'].suppress()

#Boundary conditions
def Displacement_BoundaryConditions(model, matrixName):

    mdb.models[model].DisplacementBC(amplitude=UNSET, createStepName='Initial', 
        distributionType=UNIFORM, fieldName='', localCsys=None, name='BC_X', 
        region=Region(faces=mdb.models[model].rootAssembly.instances[matrixName + '-1'].faces.getSequenceFromMask(mask=('[#1 ]', ), )), 
        u1=SET, u2=UNSET, u3=UNSET, ur1=UNSET, ur2=UNSET, ur3=UNSET)

    mdb.models[model].DisplacementBC(amplitude=UNSET, createStepName='Initial', 
        distributionType=UNIFORM, fieldName='', localCsys=None, name='BC_Y', 
        region=Region(faces=mdb.models[model].rootAssembly.instances[matrixName + '-1'].faces.getSequenceFromMask(mask=('[#8 ]', ), )), 
        u1=UNSET, u2=SET, u3=UNSET, ur1=UNSET, ur2=UNSET, ur3=UNSET)

    mdb.models[model].DisplacementBC(amplitude=UNSET, createStepName='Initial', 
        distributionType=UNIFORM, fieldName='', localCsys=None, name='BC_Z', 
        region=Region(faces=mdb.models[model].rootAssembly.instances[matrixName + '-1'].faces.getSequenceFromMask(mask=('[#20 ]', ), )), 
        u1=UNSET, u2=UNSET, u3=SET, ur1=UNSET, ur2=UNSET, ur3=UNSET)

#Temperature boundary conditions
def Temperature_BoundaryConditions(model, matrixName, temp, numGS):
    # for j in range(numGS):
    mdb.models[model].TemperatureBC(createStepName='Initial', distributionType=
        UNIFORM, fieldName='', magnitude=0.0, name='TEMP', region=Region(
        cells=
        mdb.models[model].rootAssembly.instances[matrixName + '-1'].cells.getSequenceFromMask(mask=('[#1 ]', ), )
            # +\
            # mdb.models[model].rootAssembly.instances['GS-%d-1' %(j)].cells.getSequenceFromMask(mask=('[#1 ]', ), )
        ))
    
    mdb.models[modelName].boundaryConditions['TEMP'].setValuesInStep(magnitude=temp, stepName='Step_name')

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
        initialInc=0.1, name='Step-1', noStop=OFF, 
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
    
#This function generates the mesh for the matrix
def generate_matrix_mesh(model, str_matrix, matrixMeshSize):
    
    #Mesh the matrix
    #deviationFactor and minSizeFactor have the default values from Abaqus
    mdb.models[model].parts[str_matrix].seedPart(deviationFactor=0.1, minSizeFactor=0.1, size=matrixMeshSize)
    mdb.models[model].parts[str_matrix].generateMesh()
    
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

def Create_Embedded_Elements_CuttedGS(model, matrixName, i):
        
    #Get the instance name
    str_inst = string_instance('GS', i)
    
    #Add the embedded element constraint
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

#NEED TO CREATE ONE SET PER GS VERTEX
def Create_NodeSet_Gs(model, i):

    GS_geometry = data[i]

    length_x = GS_geometry[0]
    length_y = GS_geometry[1]
    thickness = GS_geometry[2]

    part_GS = mdb.models[model].parts['GS-%d' %(i)]

    v = part_GS.vertices.getClosest(coordinates=((length_x/2,length_y/2,thickness),(-length_x/2,length_y/2,thickness),(-length_x/2,-length_y/2,thickness),
    (length_x/2,-length_y/2,thickness),(length_x/2,length_y/2,0),(-length_x/2,length_y/2,0),(-length_x/2,-length_y/2,0),(length_x/2,-length_y/2,0),)) 

    v1 = part_GS.vertices.findAt((((v[0][1])),))
    v2 = part_GS.vertices.findAt((((v[1][1])),))
    v3 = part_GS.vertices.findAt((((v[2][1])),))
    v4 = part_GS.vertices.findAt((((v[3][1])),))
    v5 = part_GS.vertices.findAt((((v[4][1])),))
    v6 = part_GS.vertices.findAt((((v[5][1])),))
    v7 = part_GS.vertices.findAt((((v[6][1])),))
    v8 = part_GS.vertices.findAt((((v[7][1])),))

    part_GS.Set(name='Vertices_GS-%d' %(i), vertices=(v1,v2,v3,v4,v5,v6,v7,v8,))
    
    
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
    #    nodes=(mdb.models['Model-1'].parts[str_matrix].nodes.getByBoundingSphere(center=P0, radius=0.001), ),
    #    name='Matrix0')
    
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
    
    #Create set
 # mdb.models['Model-1'].rootAssembly.Set(
 #     nodes=mdb.models['Model-1'].rootAssembly.instances[str_matrix+"-1"].nodes.getByBoundingSphere(center=(P0[0], P0[1]+Lxyz[1], P0[2]+Lxyz[2]), radius=0.001),
 #     name='Matrix2')
        
    #Perform union of two sets
    #mdb.models['Model-1'].rootAssembly.SetByBoolean(
    #    name='Matrix2', 
    #    operation=UNION, 
    #    sets=(mdb.models['Model-1'].rootAssembly.sets['Matrix0'], mdb.models['Model-1'].rootAssembly.sets['Matrix1'])
    #    )
    
    
#Matrix name
matrixName = 'MATRIX'
#Cutting box name
biggerBoxName = 'BIGGER_BOX'
#Name of the analysis model (unchangeable)
modelName = 'Model-1'

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

#########################################---OPENING CSV FILES---#########################################

with open(gs_file) as f:
    #x = csv.reader(f, quoting=csv.QUOTE_NONNUMERIC)
    gs_data = list(csv.reader(f, quoting=csv.QUOTE_NONNUMERIC))

N_GSs = len(gs_data)

#Read the sample geometry
with open(sample_file) as f:
    (P0, Lxyz) = list(csv.reader(f, quoting=csv.QUOTE_NONNUMERIC))

print('The number of graphene sheets is:',N_GSs)

#######################################---MODEL CREATION---#########################################


###----------------------------------Start: Create Parts----------------------------------###
#Generate the matrix part
matrix_part(P0, Lxyz, str_matrix, sheetSize)

#Box used for cutting GSs
#Create_BiggerBox(modelName, sheetSize, biggerBoxName, matrixLength)

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

#Create material and section for matrix (only elastic properties)
materials_and_sections_matrix(modelName, str_matrix, matrixMaterial, matrixSection, matrixDensity, matrixModulus, matrixPoissonR)

#Create material and section for GS (only elastic properties)
materials_and_sections_gs(N_GSs, gsMaterial, gsSection, gsDensity, gsModulus, gsPoissonR)

###----------------------------------End: Materials and sections----------------------------------###


###----------------------------------Start: Create assembly----------------------------------###
mdb.models[modelName].rootAssembly.DatumCsysByDefault(CARTESIAN)

# mdb.models[modelName].rootAssembly.Instance(
#     dependent=ON, 
#     name=biggerBoxName + '-1', 
#     part=mdb.models[modelName].parts[biggerBoxName])

#Create the cutting box using one particular GS as reference
#Create_CuttingBox(modelName, str_matrix, biggerBoxName, gs_data[1])

#Generate assembly
generate_assembly(N_GSs, modelName, str_matrix, gs_data)
###-----------------------------------End: Create assembly-----------------------------------###

# print(indexOutside)    
# print(indexInside)

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

#Create all embedded elements constraints when using option all_in
create_embedded_elments_gs(N_GSs, modelName, str_matrix, 'GS')

#Create Step and add boundary conditions
create_step_and_pbcs(modelName, P0, Lxyz, str_matrix, disp)

###-------------------------------------End: BCs-------------------------------------###


###-------------------------------------Start: Meshing-------------------------------------###

#Mesh matrix
generate_matrix_mesh(modelName, str_matrix, matrixMeshSize)

#Mesh all GS when using option all_in
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


###-------------------------------------Start: Job-------------------------------------###
#Create and submit job using Abaqus default values
mdb.Job(atTime=None, contactPrint=OFF, description='', echoPrint=OFF, 
    explicitPrecision=SINGLE, getMemoryFromAnalysis=True, historyPrint=OFF, 
    memory=90, memoryUnits=PERCENTAGE, model='Model-1', modelPrint=OFF, 
    multiprocessingMode=DEFAULT, name='Job-1', nodalOutputPrecision=SINGLE, 
    numCpus=1, numGPUs=0, queue=None, resultsFormat=ODB, scratch='', type=
    ANALYSIS, userSubroutine='', waitHours=0, waitMinutes=0)
mdb.jobs['Job-1'].submit(consistencyChecking=OFF)
###-------------------------------------End: Job--------------------------------------###

