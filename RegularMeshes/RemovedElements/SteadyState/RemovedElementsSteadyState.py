#import pdb; pdb.set_trace()

# Add Python bindings directory to PATH
import sys, os, numpy

# Initialise OpenCMISS-Iron
from opencmiss.iron import iron

# Set problem parameters
height = 1.0
width = 2.0
length = 3.0

(coordinateSystemUserNumber,
 regionUserNumber,
 basisUserNumber,
 generatedMeshUserNumber,
 mesh1UserNumber,
 mesh2UserNumber,
 decompositionUserNumber,
 geometricField1UserNumber,
 geometricField2UserNumber,
 equationsSetFieldUserNumber,
 dependentFieldUserNumber,
 equationsSetUserNumber,
 problemUserNumber) = range(1,14)

numberGlobalXElements = 5
numberGlobalYElements = 5
numberGlobalZElements = 5

iron.DiagnosticsSetOn(iron.DiagnosticTypes.IN,[1,2,3,4,5],"Diagnostics",["DOMAIN_MAPPINGS_LOCAL_FROM_GLOBAL_CALCULATE"])

numberOfComputationalNodes = iron.ComputationalNumberOfNodesGet()
computationalNodeNumber = iron.ComputationalNodeNumberGet()

mesh_component_number = 1
number_of_dimensions = 3
nodes_per_elem = 4 # for a tet mesh

# Create a RC coordinate system
coordinateSystem = iron.CoordinateSystem()
coordinateSystem.CreateStart(coordinateSystemUserNumber)
coordinateSystem.dimension = 3
coordinateSystem.CreateFinish()

# Create a region
region = iron.Region()
region.CreateStart(regionUserNumber,iron.WorldRegion)
region.label = "LaplaceRegion"
region.coordinateSystem = coordinateSystem
region.CreateFinish()

# Create a tri-linear simplex basis
basis = iron.Basis()
basis.CreateStart(basisUserNumber)
basis.TypeSet(iron.BasisTypes.SIMPLEX)
basis.numberOfXi = 3
basis.interpolationXi = [iron.BasisInterpolationSpecifications.LINEAR_SIMPLEX]*3
basis.CreateFinish()

# Create a generated mesh
generatedMesh = iron.GeneratedMesh()
generatedMesh.CreateStart(generatedMeshUserNumber,region)
generatedMesh.type = iron.GeneratedMeshTypes.REGULAR
generatedMesh.basis = [basis]
generatedMesh.extent = [width,height,length]
generatedMesh.numberOfElements = [numberGlobalXElements,numberGlobalYElements,numberGlobalZElements]

mesh1 = iron.Mesh()
generatedMesh.CreateFinish(mesh1UserNumber,mesh1)

# Find the total number of nodes and elements from the generated mesh
meshNodes = iron.MeshNodes()
mesh1.NodesGet(mesh_component_number,meshNodes)
total_number_of_nodes = meshNodes.NumberOfNodesGet()

total_number_of_elements = mesh1.NumberOfElementsGet()

# Start the creation of a second manually generated mesh
mesh2 = iron.Mesh()
mesh2.CreateStart(mesh2UserNumber,region,number_of_dimensions)
mesh2.NumberOfComponentsSet(1)

# Define which elements are included in the new mesh, and the total number of elements
original_elements_list = range(1,total_number_of_elements+1)
removed_elements = (numpy.array([1,750]) - 1).tolist() # list of elements to be removed from mesh
new_elements_list = numpy.delete(original_elements_list, removed_elements).tolist()

mesh2.NumberOfElementsSet(len(new_elements_list))

# Define elements from the original mesh
originalElements = iron.MeshElements()
mesh1.ElementsGet(mesh_component_number, originalElements)

# Find nodes associated with a particular element from original mesh
elements = iron.MeshElements()
elements.CreateStart(mesh2,mesh_component_number,basis)
selected_element_nodes = []

# Copy the desired nodal values to second mesh
for idx, elem_num in enumerate(new_elements_list):
    element_nodes = originalElements.NodesGet(elem_num,nodes_per_elem)
    elements.NodesSet(idx+1,element_nodes)
    selected_element_nodes += element_nodes.tolist()

# Refers to elements by their user number as described in the original mesh
elements.UserNumbersAllSet(new_elements_list)

# Removes duplicates from the list of nodal values and orders them
ordered_nodes = dict.fromkeys(selected_element_nodes).keys()

elements.CreateFinish()

mesh2.CreateFinish()

# Create a decomposition for the mesh
decomposition = iron.Decomposition()
decomposition.CreateStart(decompositionUserNumber,mesh2)
decomposition.type = iron.DecompositionTypes.CALCULATED
decomposition.numberOfDomains = numberOfComputationalNodes
decomposition.CreateFinish()

# Create a field for the geometry
geometricField1 = iron.Field()
geometricField1.CreateStart(geometricField1UserNumber,region)
geometricField1.meshDecomposition = decomposition
geometricField1.ComponentMeshComponentSet(iron.FieldVariableTypes.U,1,1)
geometricField1.ComponentMeshComponentSet(iron.FieldVariableTypes.U,2,1)
geometricField1.ComponentMeshComponentSet(iron.FieldVariableTypes.U,3,1)
geometricField1.CreateFinish()

generatedMesh.GeometricParametersCalculate(geometricField1)

# Create a field for the geometry
geometricField2 = iron.Field()
geometricField2.CreateStart(geometricField2UserNumber,region)
geometricField2.meshDecomposition = decomposition
geometricField2.ComponentMeshComponentSet(iron.FieldVariableTypes.U,1,1)
geometricField2.ComponentMeshComponentSet(iron.FieldVariableTypes.U,2,1)
geometricField2.ComponentMeshComponentSet(iron.FieldVariableTypes.U,3,1)
geometricField2.CreateFinish()

# Update the geometric field parameters
geometricField2.ParameterSetUpdateStart(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)

#for node_num in range (1,len(ordered_nodes)+1):
for node_num in ordered_nodes:
    field_val_1 = geometricField1.ParameterSetGetNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,1,node_num,1)
    field_val_2 = geometricField1.ParameterSetGetNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,1,node_num,2)
    field_val_3 = geometricField1.ParameterSetGetNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,1,node_num,3)

    geometricField2.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,1,node_num,1,field_val_1)
    geometricField2.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,1,node_num,2,field_val_2)
    geometricField2.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,1,node_num,3,field_val_3)

geometricField2.ParameterSetUpdateFinish(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)

# Create standard Laplace equations set
equationsSetField = iron.Field()
equationsSet = iron.EquationsSet()
equationsSetSpecification = [iron.EquationsSetClasses.CLASSICAL_FIELD,
                             iron.EquationsSetTypes.LAPLACE_EQUATION,
                             iron.EquationsSetSubtypes.STANDARD_LAPLACE]
equationsSet.CreateStart(equationsSetUserNumber,region,geometricField2,
                         equationsSetSpecification,equationsSetFieldUserNumber,equationsSetField)
equationsSet.CreateFinish()

# Create dependent field
dependentField = iron.Field()
equationsSet.DependentCreateStart(dependentFieldUserNumber,dependentField)
dependentField.DOFOrderTypeSet(iron.FieldVariableTypes.U,iron.FieldDOFOrderTypes.SEPARATED)
dependentField.DOFOrderTypeSet(iron.FieldVariableTypes.DELUDELN,iron.FieldDOFOrderTypes.SEPARATED)
equationsSet.DependentCreateFinish()

# Initialise dependent field
dependentField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,0.5)

# Create equations
equations = iron.Equations()
equationsSet.EquationsCreateStart(equations)
equations.sparsityType = iron.EquationsSparsityTypes.SPARSE
equations.outputType = iron.EquationsOutputTypes.NONE
equationsSet.EquationsCreateFinish()

# Create problem
problem = iron.Problem()
problemSpecification = [iron.ProblemClasses.CLASSICAL_FIELD,
                        iron.ProblemTypes.LAPLACE_EQUATION,
                        iron.ProblemSubtypes.STANDARD_LAPLACE]
problem.CreateStart(problemUserNumber, problemSpecification)
problem.CreateFinish()

# Create control loops
problem.ControlLoopCreateStart()
problem.ControlLoopCreateFinish()

# Create problem solver
solver = iron.Solver()
problem.SolversCreateStart()
problem.SolverGet([iron.ControlLoopIdentifiers.NODE],1,solver)
solver.outputType = iron.SolverOutputTypes.SOLVER
solver.linearType = iron.LinearSolverTypes.ITERATIVE
solver.linearIterativeAbsoluteTolerance = 1.0E-12
solver.linearIterativeRelativeTolerance = 1.0E-12
problem.SolversCreateFinish()

# Create solver equations and add equations set to solver equations
solver = iron.Solver()
solverEquations = iron.SolverEquations()
problem.SolverEquationsCreateStart()
problem.SolverGet([iron.ControlLoopIdentifiers.NODE],1,solver)
solver.SolverEquationsGet(solverEquations)
solverEquations.sparsityType = iron.SolverEquationsSparsityTypes.SPARSE
equationsSetIndex = solverEquations.EquationsSetAdd(equationsSet)
problem.SolverEquationsCreateFinish()

# Create boundary conditions and set first and last nodes to 0.0 and 1.0
boundaryConditions = iron.BoundaryConditions()
solverEquations.BoundaryConditionsCreateStart(boundaryConditions)
firstNodeNumber=1
nodes = iron.Nodes()
region.NodesGet(nodes)
lastNodeNumber = nodes.numberOfNodes
firstNodeDomain = decomposition.NodeDomainGet(firstNodeNumber,1)
lastNodeDomain = decomposition.NodeDomainGet(lastNodeNumber,1)

if firstNodeDomain == computationalNodeNumber:
    boundaryConditions.SetNode(dependentField,iron.FieldVariableTypes.U,1,1,firstNodeNumber,1,iron.BoundaryConditionsTypes.FIXED,0.0)
if lastNodeDomain == computationalNodeNumber:
    boundaryConditions.SetNode(dependentField,iron.FieldVariableTypes.U,1,1,lastNodeNumber,1,iron.BoundaryConditionsTypes.FIXED,1.0)
solverEquations.BoundaryConditionsCreateFinish()

# Solve the problem
problem.Solve()

# Export results as fml files for the geometric field and the dependent field (.geometric and .phi respectively). Outputs a .xml file.
baseName = "laplace"
dataFormat = "PLAIN_TEXT"
fml = iron.FieldMLIO()
fml.OutputCreate(mesh2, "", baseName, dataFormat)
fml.OutputAddFieldNoType(baseName+".geometric", dataFormat, geometricField2,
                         iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES)
fml.OutputAddFieldNoType(baseName+".phi", dataFormat, dependentField,
                         iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES)
fml.OutputWrite("RemovedElementsSteadyState.xml")
fml.Finalise()

# Export results
fields = iron.Fields()
fields.CreateRegion(region)
fields.NodesExport("RemovedElementsSteadyStateResults","FORTRAN")
fields.ElementsExport("RemovedElementsSteadyStateResults","FORTRAN")
fields.Finalise()

iron.Finalise()