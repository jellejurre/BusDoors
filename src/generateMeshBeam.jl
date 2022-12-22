using Gmsh
import Gmsh: gmsh
function MeshGeneratorBeam(x1, y1, z1, meshSize)
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 1)
    gmsh.option.setNumber("Mesh.Algorithm", 6)
    gmsh.clear()
    gmsh.model.add("geometry")

    # Add points
    gmsh.model.geo.addPoint(0,   0,   0, meshSize, 11)
    gmsh.model.geo.addPoint(0,   y1,   0, meshSize, 12)
    gmsh.model.geo.addPoint(x1,   y1,   0, meshSize, 13)
    gmsh.model.geo.addPoint(x1,   0,   0, meshSize, 14)
    gmsh.model.geo.addPoint(0,   0,   z1, meshSize, 15)
    gmsh.model.geo.addPoint(0,   y1,   z1, meshSize, 16)
    gmsh.model.geo.addPoint(x1,   y1,   z1, meshSize, 17)
    gmsh.model.geo.addPoint(x1,   0,   z1, meshSize, 18)
    
    # Add lines
    gmsh.model.geo.addLine( 11,  12,  11)
    gmsh.model.geo.addLine( 12,  13,  12)
    gmsh.model.geo.addLine( 13,  14,  13)
    gmsh.model.geo.addLine( 11,  14,  14)
    gmsh.model.geo.addLine( 15,  16,  15)
    gmsh.model.geo.addLine( 16,  17,  16)
    gmsh.model.geo.addLine( 17,  18,  17)
    gmsh.model.geo.addLine( 15,  18,  18)
    gmsh.model.geo.addLine( 15,  11,  19)
    gmsh.model.geo.addLine( 16,  12,  110)
    gmsh.model.geo.addLine( 17,  13,  111)
    gmsh.model.geo.addLine( 18,  14,  112)

    # Construct curve loops and surfaces    (watch out for line direction)
    gmsh.model.geo.addCurveLoop([11, 12, 13, -14], 11)
    gmsh.model.geo.addCurveLoop([15, 16, 17, -18], 12)
    gmsh.model.geo.addCurveLoop([11, -110, -15, 19], 13)
    gmsh.model.geo.addCurveLoop([-13, -111, 17, 112], 14)
    gmsh.model.geo.addCurveLoop([-14, -19, 18, 112], 15)
    gmsh.model.geo.addCurveLoop([111, -12, -110, 16], 16)

    gmsh.model.geo.addPlaneSurface([11], 11)
    gmsh.model.geo.addPlaneSurface([12], 12)
    gmsh.model.geo.addPlaneSurface([13], 13)
    gmsh.model.geo.addPlaneSurface([14], 14)
    gmsh.model.geo.addPlaneSurface([15], 15)
    gmsh.model.geo.addPlaneSurface([16], 16)

    # outer volume aluminium
    gmsh.model.geo.addSurfaceLoop([11, 12, 13, 14, 15, 16], 11)
    gmsh.model.geo.addVolume([11], 11)

    # Physical groups
    # def: addPhysicalGroup(dim, tags, tag = -1, name = "")
    #   Add a physical group of dimension `dim`, grouping the model entities with tags
    #   `tags`. Return the tag of the physical group, equal to `tag` if `tag` is
    #   positive, or a new tag if `tag` < 0. Set the name of the physical group if
    #   `name` is not empty.

    gmsh.model.addPhysicalGroup(1, [11,12,13,14,15, 16, 17, 18, 19, 110, 111, 112], 11)
    gmsh.model.setPhysicalName(1, 11, "FreeEdges1")
    gmsh.model.addPhysicalGroup(2, [11, 12, 13, 15, 16], 21)
    gmsh.model.setPhysicalName(2, 21, "FreeAreas1")
    gmsh.model.addPhysicalGroup(2, [14], 22)
    gmsh.model.setPhysicalName(2, 22, "Dirichlet1")
    gmsh.model.addPhysicalGroup(3, [11], 31)
    gmsh.model.setPhysicalName(3, 31, "Beam1")


    gmsh.model.geo.synchronize()
    # We can then generate a 3D mesh...
    gmsh.model.mesh.generate(3)
    # ... and save it to disk
    gmsh.write("geometry_beam.msh")
    gmsh.finalize()
end

x1 = 1
y1 = 0.05
z1 = 0.05
MeshGeneratorBeam(x1, y1, z1, 0.01)