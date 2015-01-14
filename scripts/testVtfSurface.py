from GoTools import *

# First block - unit square
surf1 = Surface(2, 2, [0,0,1,1], [0,0,1,1],
                [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [1.0, 1.0, 0.0]])

# Second block - translated unit square
surf2 = Surface(2, 2, [0,0,1,1], [0,0,1,1],
                [[1.0, 0.0, 0.0], [2.0, 0.0, 0.0], [1.0, 1.0, 0.0], [2.0, 1.0, 0.0]])

# First data field
surf1s = surf1.Clone([0.0, 1.0, 0.0, 1.0])

# Second data field
surf2s = surf2.Clone([1.0, 2.0, 1.0, 2.0])

# Open VTF file
vtf = VTF("/tmp/regtest.g2",binary=False)

# Tesselate first surface
(nod, elem) = surf1.Tesselate(8,8)
# Write FEM model to VTF file
vtf.AddGeometryBlock(nod, elem, 1, 2)

# Tesselate second surface
(nod2, elem2) = surf2.Tesselate()
# Write FEM model to VTF file
vtf.AddGeometryBlock(nod2, elem2, 2, 2)

# Write out grid descriptor
vtf.AddGeometryDescriptor(2)

# Evaluate fields
(t11, t21) = surf1s.GetTesselationParams(8,8)
data1 = surf1s.EvaluateGrid(t11, t21)
(t12, t22) = surf2s.GetTesselationParams()
data2 = surf2s.EvaluateGrid(t12, t22)

# Write fields to file
femfields = [vtf.AddField(data1, 1), vtf.AddField(data2, 2)]
vtf.AddFieldBlocks(femfields, "test")

# Displacement fields
surf1d = surf1.Clone([[1.0, 0.0], [1.0, 0.0], [1.0, 0.0], [1.0, 0.0]])
data1d = surf1d.EvaluateGrid(t11, t21)
surf2d = surf2.Clone([[1.0, 0.0], [1.0, 0.0], [1.0, 0.0], [1.0, 0.0]])
data2d = surf2d.EvaluateGrid(t12, t22)

# Write displacement fieds to file
femfields_d = [vtf.AddField(data1d, 1), vtf.AddField(data2d, 2)]
vtf.AddFieldBlocks(femfields_d, "test_v", 2, True)

# Vector fields
surf1v = surf1.Clone([[0.0, 1.0], [0.0, 1.0], [0.0, 1.0], [0.0, 1.0]])
data1v = surf1v.EvaluateGrid(t11, t21)
surf2v = surf2.Clone([[1.0, 0.0], [1.0, 0.0], [1.0, 0.0], [1.0, 0.0]])
data2v = surf2v.EvaluateGrid(t12, t22)

femfields_v = [vtf.AddField(data1v, 1), vtf.AddField(data2v, 2)]

vtf.AddFieldBlocks(femfields_v, "test_v2", 2)
vtf.AddState(0.0)

# Second time step
surf1s = surf1.Clone([0.0, 5.0, 0.0, 5.0])
surf2s = surf2.Clone([2.0, 3.0, 2.0, 3.0])
data1 = surf1s.EvaluateGrid(t11, t21)
data2 = surf2s.EvaluateGrid(t12, t22)

femfields = [vtf.AddField(data1, 1), vtf.AddField(data2, 2)]
vtf.AddFieldBlocks(femfields, "test")
vtf.AddState(1.0)
