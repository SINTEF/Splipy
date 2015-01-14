from GoTools import *

# First block - unit square
crv1 = Curve(2, [0,0,1,1], [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]])

# Second block - translated unit square
crv2 = Curve(2, [0,0,1,1], [[1.0, 0.0, 0.0],[2.0, 0.0, 0.0]])

# First data field
crv1s = crv1.Clone([0.0, 1.0])

# Second data field
crv2s = crv2.Clone([1.0, 2.0])

# Open VTF file
vtf = VTF("/tmp/regtest.g2",binary=False)

# Tesselate first curve
(nod, elem) = crv1.Tesselate(8)
# Write FEM model to VTF file
vtf.AddGeometryBlock(nod, elem, 1, 1)

# Tesselate second curve
(nod2, elem2) = crv2.Tesselate()
# Write FEM model to VTF file
vtf.AddGeometryBlock(nod2, elem2, 2, 1)

# Write out grid descriptor
vtf.AddGeometryDescriptor(2)

# Evaluate fields
t1 = crv1s.GetTesselationParams(8)
data1 = crv1s.EvaluateGrid(t1)
t2 = crv2s.GetTesselationParams()
data2 = crv2s.EvaluateGrid(t2)

# Write fields to file
femfields = [vtf.AddField(data1, 1), vtf.AddField(data2, 2)]
vtf.AddFieldBlocks(femfields, "test")

vtf.AddState(0.0)

# Second time step
crv1s = crv1.Clone([0.0, 5.0])
crv2s = crv2.Clone([2.0, 3.0])
data1 = crv1s.EvaluateGrid(t1)
data2 = crv2s.EvaluateGrid(t2)

femfields = [vtf.AddField(data1, 1), vtf.AddField(data2, 2)]
vtf.AddFieldBlocks(femfields, "test")
vtf.AddState(1.0)
