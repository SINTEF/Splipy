from GoTools import *

# First block - unit cube
vol1 = Volume(2, 2, 2, [0,0,1,1], [0,0,1,1], [0,0,1,1],
               [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [1.0, 1.0, 0.0],
                [0.0, 0.0, 1.0], [1.0, 0.0, 1.0], [0.0, 1.0, 1.0], [1.0, 1.0, 1.0]])

# Second block - translated unit cube
vol2  = Volume(2, 2, 2, [0,0,1,1], [0,0,1,1], [0,0,1,1],
               [[1.0, 0.0, 0.0], [2.0, 0.0, 0.0], [1.0, 1.0, 0.0], [2.0, 1.0, 0.0],
                [1.0, 0.0, 1.0], [2.0, 0.0, 1.0], [1.0, 1.0, 1.0], [2.0, 1.0, 1.0]])

# First data field
vol1s = vol1.Clone([0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0])

# Second data field
vol2s = vol2.Clone([1.0, 2.0, 1.0, 2.0, 1.0, 2.0, 1.0, 2.0])

# Open VTF file
vtf = VTF("/tmp/regtest.g2",binary=False)

# Tesselate first surface
(nod, elem) = vol1.Tesselate(8,8,8)
# Write FEM model to VTF file
vtf.AddGeometryBlock(nod, elem, 1, 3)

# Tesselate second surface
(nod2, elem2) = vol2.Tesselate()
# Write FEM model to VTF file
vtf.AddGeometryBlock(nod2, elem2, 2, 3)

# Write out grid descriptor
vtf.AddGeometryDescriptor(2)

# Evaluate fields
(t11, t21, t31) = vol1s.GetTesselationParams(8,8,8)
data1 = vol1s.EvaluateGrid(t11, t21, t31)
(t12, t22, t32) = vol2s.GetTesselationParams()
data2 = vol2s.EvaluateGrid(t12, t22, t32)

# Write fields to file
femfields = [vtf.AddField(data1, 1), vtf.AddField(data2, 2)]
vtf.AddFieldBlocks(femfields, "test")

vtf.AddState(0.0)
