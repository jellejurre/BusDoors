using Gmsh
import Gmsh: gmsh
function MeshGeneratorBeam(x1, y, z1, x2, z2, meshSize)
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 1)
    gmsh.option.setNumber("Mesh.Algorithm", 6)
    gmsh.clear()
    gmsh.model.add("geometry")

    # Add points
    gmsh.model.geo.addPoint(x1,   0,   z1, meshSize, 1)
    gmsh.model.geo.addPoint(x1,   0,   0, meshSize, 2)
    gmsh.model.geo.addPoint(0,   0,   0, meshSize, 3)
    gmsh.model.geo.addPoint(0,   0,   z2, meshSize, 4)
    gmsh.model.geo.addPoint(x2,   0,   z2, meshSize, 5)
    gmsh.model.geo.addPoint(x2,   0,   z1, meshSize, 6)
    gmsh.model.geo.addPoint(x1,   y,   z1, meshSize, 7)
    gmsh.model.geo.addPoint(x1,   y,   0, meshSize, 8)
    gmsh.model.geo.addPoint(0,   y,   0, meshSize, 9)
    gmsh.model.geo.addPoint(0,   y,   z2, meshSize, 10)
    gmsh.model.geo.addPoint(x2,   y,   z2, meshSize, 11)
    gmsh.model.geo.addPoint(x2,   y,   z1, meshSize, 12)
    
    # Add lines
    gmsh.model.geo.addLine(1, 2, 1)
    gmsh.model.geo.addLine(2, 3, 2)
    gmsh.model.geo.addLine(3, 4, 3)
    gmsh.model.geo.addLine(4, 5, 4)
    gmsh.model.geo.addLine(5, 6, 5)
    gmsh.model.geo.addLine(1, 6, 6)
    gmsh.model.geo.addLine(2, 8, 7)
    gmsh.model.geo.addLine(3, 9, 8)
    gmsh.model.geo.addLine(4, 10, 9)
    gmsh.model.geo.addLine(5, 11, 10)
    gmsh.model.geo.addLine(6, 12, 11)
    gmsh.model.geo.addLine(1, 7, 12)
    gmsh.model.geo.addLine(7, 8, 13)
    gmsh.model.geo.addLine(8, 9, 14)
    gmsh.model.geo.addLine(9, 10, 15)
    gmsh.model.geo.addLine(10, 11, 16)
    gmsh.model.geo.addLine(11, 12, 17)
    gmsh.model.geo.addLine(7, 12, 18)

    # Construct curve loops and surfaces    (watch out for line direction)
    gmsh.model.geo.addCurveLoop([1, 2, 3, 4, 5, -6], 101)
    gmsh.model.geo.addCurveLoop([1, 7, -13, -12], 102)
    gmsh.model.geo.addCurveLoop([2, 8, -14, -7], 103)
    gmsh.model.geo.addCurveLoop([3, 9, -15, -8], 104)
    gmsh.model.geo.addCurveLoop([4, 10, -16, -9], 105)
    gmsh.model.geo.addCurveLoop([5, 11, -17, -10], 106)
    gmsh.model.geo.addCurveLoop([-6, 12, 18, -11], 107)
    gmsh.model.geo.addCurveLoop([13, 14, 15, 16, 17, -18], 108)

    gmsh.model.geo.addPlaneSurface([101], 201)
    gmsh.model.geo.addPlaneSurface([102], 202)
    gmsh.model.geo.addPlaneSurface([103], 203)
    gmsh.model.geo.addPlaneSurface([104], 204)
    gmsh.model.geo.addPlaneSurface([105], 205)
    gmsh.model.geo.addPlaneSurface([106], 206)
    gmsh.model.geo.addPlaneSurface([107], 207)
    gmsh.model.geo.addPlaneSurface([108], 208)

    # outer volume aluminium
    gmsh.model.geo.addSurfaceLoop([201, 202, 203, 204, 205, 206, 207, 208], 301)
    gmsh.model.geo.addVolume([301], 401)

    # Physical groups
    # def: addPhysicalGroup(dim, tags, tag = -1, name = "")
    #   Add a physical group of dimension `dim`, grouping the model entities with tags
    #   `tags`. Return the tag of the physical group, equal to `tag` if `tag` is
    #   positive, or a new tag if `tag` < 0. Set the name of the physical group if
    #   `name` is not empty.

    gmsh.model.addPhysicalGroup(1, [1,2,3,4,5,6,7,8,9,1011,12,13,14,15, 16, 17, 18], 11)
    gmsh.model.setPhysicalName(1, 11, "FreeEdges")
    gmsh.model.addPhysicalGroup(2, [201, 203, 204, 206, 207, 208], 21)
    gmsh.model.setPhysicalName(2, 21, "FreeAreas")
    gmsh.model.addPhysicalGroup(2, [202], 22)
    gmsh.model.setPhysicalName(2, 22, "BeamEnd1")
    gmsh.model.addPhysicalGroup(2, [205], 23)
    gmsh.model.setPhysicalName(2, 23, "BeamEnd2")
    gmsh.model.addPhysicalGroup(3, [401], 31)
    gmsh.model.setPhysicalName(3, 31, "Beam")


    gmsh.model.geo.synchronize()
    # We can then generate a 3D mesh...
    gmsh.model.mesh.generate(3)
    # ... and save it to disk
    gmsh.write("geometry_beam2.msh")
    gmsh.finalize()
end

x1 = 1
y1 = 0.05
z1 = 0.05
x2 = 0.05
z2 = 1
MeshGeneratorBeam(x1, y1, z1, x2, z2, 0.01)