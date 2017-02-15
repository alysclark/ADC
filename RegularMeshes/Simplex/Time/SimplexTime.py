# import pdb; pdb.set_trace()

# Add Python bindings directory to PATH
import sys, os, numpy

# Initialise OpenCMISS-Iron
from opencmiss.iron import iron

# Set problem parameters
height = 0.1
width = 0.2
length = 0.3
diff_coeff = 22.5  # in mm^2/sec
initial_conc = 0.5
start_time = 0.0
end_time = 10000
time_step = 10
screen_output_freq = 2  # how many time steps between outputs to screen

(coordinateSystemUserNumber,
 regionUserNumber,
 basisUserNumber,
 generatedMeshUserNumber,
 meshUserNumber,
 decompositionUserNumber,
 geometricFieldUserNumber,
 equationsSetFieldUserNumber,
 dependentFieldUserNumber,
 materialFieldUserNumber,
 equationsSetUserNumber,
 problemUserNumber) = range(1, 13)

numberGlobalXElements = 5
numberGlobalYElements = 5
numberGlobalZElements = 5

iron.DiagnosticsSetOn(iron.DiagnosticTypes.IN, [1, 2, 3, 4, 5], "Diagnostics",
                      ["DOMAIN_MAPPINGS_LOCAL_FROM_GLOBAL_CALCULATE"])

# Get the computational nodes information
numberOfComputationalNodes = iron.ComputationalNumberOfNodesGet()
computationalNodeNumber = iron.ComputationalNodeNumberGet()

# Create a RC coordinate system
coordinateSystem = iron.CoordinateSystem()
coordinateSystem.CreateStart(coordinateSystemUserNumber)
coordinateSystem.dimension = 3
coordinateSystem.CreateFinish()

# Create a region
region = iron.Region()
region.CreateStart(regionUserNumber, iron.WorldRegion)
region.label = "LaplaceRegion"
region.coordinateSystem = coordinateSystem
region.CreateFinish()

# Create a tri-linear simplex basis
basis = iron.Basis()
basis.CreateStart(basisUserNumber)
basis.TypeSet(iron.BasisTypes.SIMPLEX)
basis.numberOfXi = 3
basis.interpolationXi = [iron.BasisInterpolationSpecifications.LINEAR_SIMPLEX] * 3
basis.CreateFinish()

# Create a generated mesh
generatedMesh = iron.GeneratedMesh()
generatedMesh.CreateStart(generatedMeshUserNumber, region)
generatedMesh.type = iron.GeneratedMeshTypes.REGULAR
generatedMesh.basis = [basis]
generatedMesh.extent = [width, height, length]
generatedMesh.numberOfElements = [numberGlobalXElements, numberGlobalYElements, numberGlobalZElements]

mesh = iron.Mesh()
generatedMesh.CreateFinish(meshUserNumber, mesh)

# Create a decomposition for the mesh
decomposition = iron.Decomposition()
decomposition.CreateStart(decompositionUserNumber, mesh)
decomposition.type = iron.DecompositionTypes.CALCULATED
decomposition.numberOfDomains = numberOfComputationalNodes
decomposition.CreateFinish()

# Create a field for the geometry
geometricField = iron.Field()
geometricField.CreateStart(geometricFieldUserNumber, region)
geometricField.meshDecomposition = decomposition
geometricField.ComponentMeshComponentSet(iron.FieldVariableTypes.U, 1, 1)
geometricField.ComponentMeshComponentSet(iron.FieldVariableTypes.U, 2, 1)
geometricField.ComponentMeshComponentSet(iron.FieldVariableTypes.U, 3, 1)
geometricField.CreateFinish()

# Set geometry from the generated mesh
generatedMesh.GeometricParametersCalculate(geometricField)

# Create standard Laplace equations set
equationsSetField = iron.Field()
equationsSet = iron.EquationsSet()
equationsSetSpecification = [iron.EquationsSetClasses.CLASSICAL_FIELD,
                             iron.EquationsSetTypes.DIFFUSION_EQUATION,
                             iron.EquationsSetSubtypes.NO_SOURCE_DIFFUSION]
equationsSet.CreateStart(equationsSetUserNumber, region, geometricField,
                         equationsSetSpecification, equationsSetFieldUserNumber, equationsSetField)
equationsSet.CreateFinish()

# Create dependent field
dependentField = iron.Field()
equationsSet.DependentCreateStart(dependentFieldUserNumber, dependentField)
dependentField.DOFOrderTypeSet(iron.FieldVariableTypes.U, iron.FieldDOFOrderTypes.SEPARATED)
dependentField.DOFOrderTypeSet(iron.FieldVariableTypes.DELUDELN, iron.FieldDOFOrderTypes.SEPARATED)
equationsSet.DependentCreateFinish()

# Create material field
materialField = iron.Field()
equationsSet.MaterialsCreateStart(materialFieldUserNumber, materialField)

# Sets the material field component number
materialField.ComponentMeshComponentSet(iron.FieldVariableTypes.U, 1, 1)
materialField.ComponentMeshComponentSet(iron.FieldVariableTypes.U, 2, 1)
materialField.ComponentMeshComponentSet(iron.FieldVariableTypes.U, 3, 1)

# Change to nodal based interpolation
materialField.ComponentInterpolationSet(iron.FieldVariableTypes.U, 1, iron.FieldInterpolationTypes.NODE_BASED)
materialField.ComponentInterpolationSet(iron.FieldVariableTypes.U, 2, iron.FieldInterpolationTypes.NODE_BASED)
materialField.ComponentInterpolationSet(iron.FieldVariableTypes.U, 3, iron.FieldInterpolationTypes.NODE_BASED)

equationsSet.MaterialsCreateFinish()

# Changing diffusion coefficient
materialField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES, 1, diff_coeff)
materialField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES, 2, diff_coeff)
materialField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES, 3, diff_coeff)

# Initialise dependent field
dependentField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES, 1,
                                           initial_conc)

# Create equations
equations = iron.Equations()
equationsSet.EquationsCreateStart(equations)
equations.sparsityType = iron.EquationsSparsityTypes.SPARSE
equations.outputType = iron.EquationsOutputTypes.NONE
equationsSet.EquationsCreateFinish()

# Create problem
problem = iron.Problem()
problemSpecification = [iron.ProblemClasses.CLASSICAL_FIELD,
                        iron.ProblemTypes.DIFFUSION_EQUATION,
                        iron.ProblemSubtypes.NO_SOURCE_DIFFUSION]
problem.CreateStart(problemUserNumber, problemSpecification)
problem.CreateFinish()

# Create control loops
problem.ControlLoopCreateStart()
controlLoop = iron.ControlLoop()
problem.ControlLoopGet([iron.ControlLoopIdentifiers.NODE], controlLoop)
controlLoop.TimeOutputSet(screen_output_freq)
problem.ControlLoopCreateFinish()

# Create problem solver
dynamicSolver = iron.Solver()
problem.SolversCreateStart()
problem.SolverGet([iron.ControlLoopIdentifiers.NODE], 1, dynamicSolver)
dynamicSolver.outputType = iron.SolverOutputTypes.NONE
linearSolver = iron.Solver()
dynamicSolver.DynamicLinearSolverGet(linearSolver)
linearSolver.outputType = iron.SolverOutputTypes.NONE
linearSolver.linearType = iron.LinearSolverTypes.DIRECT
# linearSolver.LinearIterativeMaximumIterationsSet(1000)
problem.SolversCreateFinish()

# Create solver equations and add equations set to solver equations
solver = iron.Solver()
solverEquations = iron.SolverEquations()
problem.SolverEquationsCreateStart()
problem.SolverGet([iron.ControlLoopIdentifiers.NODE], 1, solver)
solver.SolverEquationsGet(solverEquations)
solverEquations.sparsityType = iron.SolverEquationsSparsityTypes.SPARSE
equationsSetIndex = solverEquations.EquationsSetAdd(equationsSet)
problem.SolverEquationsCreateFinish()

# Create boundary conditions and set first and last nodes to 0.0 and 1.0
boundaryConditions = iron.BoundaryConditions()
solverEquations.BoundaryConditionsCreateStart(boundaryConditions)
firstNodeNumber = 1
nodes = iron.Nodes()
region.NodesGet(nodes)
lastNodeNumber = nodes.numberOfNodes
firstNodeDomain = decomposition.NodeDomainGet(firstNodeNumber, 1)
lastNodeDomain = decomposition.NodeDomainGet(lastNodeNumber, 1)
if firstNodeDomain == computationalNodeNumber:
    boundaryConditions.SetNode(dependentField, iron.FieldVariableTypes.U, 1, 1, firstNodeNumber, 1,
                               iron.BoundaryConditionsTypes.FIXED, 0.0)
if lastNodeDomain == computationalNodeNumber:
    boundaryConditions.SetNode(dependentField, iron.FieldVariableTypes.U, 1, 1, lastNodeNumber, 1,
                               iron.BoundaryConditionsTypes.FIXED, 1.0)
solverEquations.BoundaryConditionsCreateFinish()

# While loop solving at different time steps
current_field_array = numpy.zeros(216)
previous_field_array = numpy.zeros(216)

# Set time-dependent parameters
number_of_steps = 0
tolerance_met = 0

while tolerance_met == 0:
    number_of_steps += 1

    # Set the new time loop
    print 'Control loop:', number_of_steps
    print 'Time step:', time_step
    print 'Start time:', start_time
    print 'End time:', end_time

    # Get the previous field values
    if number_of_steps != 1:
        for node_num in range(1, lastNodeNumber + 1):
            previous_field_array[(node_num) - 1] = dependentField.ParameterSetGetNodeDP(iron.FieldVariableTypes.U,
                                                                                       iron.FieldParameterSetTypes.PREVIOUS_VALUES,
                                                                                       1, 1,
                                                                                       node_num, 1)

    # Solve the problem using the new time loop
    problem.Solve()

    # Get the current field values
    for node_num in range(1, lastNodeNumber + 1):
        current_field_array[(node_num) - 1] = dependentField.ParameterSetGetNodeDP(iron.FieldVariableTypes.U,
                                                                                   iron.FieldParameterSetTypes.VALUES,
                                                                                   1, 1,
                                                                                   node_num, 1)

    print 'The field array for loop', number_of_steps, 'is:\n', current_field_array

    difference = max(abs(current_field_array - previous_field_array))
    # tolerance = min((1E-02 + 1E-01 * max(abs(current_field_array))), (1E-04 + 1E-04 * max(abs(previous_field_array))))
    tolerance = 1E-03

    print '\nThe maximum difference between the current and previous field is:', difference

    print '\nThe tolerance to meet is:', tolerance

    # Export results
    fields = iron.Fields()
    fields.CreateRegion(region)
    fields.NodesExport("SimplexTimeResults", "FORTRAN")
    fields.ElementsExport("SimplexTimeResults", "FORTRAN")
    fields.Finalise()

    if number_of_steps != 1:
        if difference <= tolerance:
            tolerance_met = 1

print "Total number of iterations: ", number_of_steps

iron.Finalise()