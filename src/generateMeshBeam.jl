using Gmsh
import Gmsh: gmsh
function MeshGeneratorBeam(x1, y1, z1, x2, y2, z2, second_beam, meshSize)
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

    if second_beam
        
        # low and high points of x2, y2 and z2
        x2l = 0
        x2h = x2
        y2l = y1
        y2h = y1 + y2
        z2l = 0
        z2h = z2

        gmsh.model.geo.addPoint(x2l,   y2l,   z2l, meshSize, 21)
        gmsh.model.geo.addPoint(x2l,   y2h,   z2l, meshSize, 22)
        gmsh.model.geo.addPoint(x2h,   y2h,   z2l, meshSize, 23)
        gmsh.model.geo.addPoint(x2h,   y2l,   z2l, meshSize, 24)
        gmsh.model.geo.addPoint(x2l,   y2l,   z2h, meshSize, 25)
        gmsh.model.geo.addPoint(x2l,   y2h,   z2h, meshSize, 26)
        gmsh.model.geo.addPoint(x2h,   y2h,   z2h, meshSize, 27)
        gmsh.model.geo.addPoint(x2h,   y2l,   z2h, meshSize, 28)

        gmsh.model.geo.addLine( 21,  22,  21)
        gmsh.model.geo.addLine( 22,  23,  22)
        gmsh.model.geo.addLine( 23,  24,  23)
        gmsh.model.geo.addLine( 21,  24,  24)
        gmsh.model.geo.addLine( 25,  26,  25)
        gmsh.model.geo.addLine( 26,  27,  26)
        gmsh.model.geo.addLine( 27,  28,  27)
        gmsh.model.geo.addLine( 25,  28,  28)
        gmsh.model.geo.addLine( 25,  21,  29)
        gmsh.model.geo.addLine( 26,  22,  210)
        gmsh.model.geo.addLine( 27,  23,  211)
        gmsh.model.geo.addLine( 28,  24,  212)

        gmsh.model.geo.addCurveLoop([21, 22, 23, -24], 21)
        gmsh.model.geo.addCurveLoop([25, 26, 27, -28], 22)
        gmsh.model.geo.addCurveLoop([21, -210, -25, 29], 23)
        gmsh.model.geo.addCurveLoop([-23, -211, 27, 212], 24)
        gmsh.model.geo.addCurveLoop([-24, -29, 28, 212], 25)
        gmsh.model.geo.addCurveLoop([211, -22, -210, 26], 26)

        gmsh.model.geo.addPlaneSurface([21], 21)
        gmsh.model.geo.addPlaneSurface([22], 22)    
        gmsh.model.geo.addPlaneSurface([23], 23)    
        gmsh.model.geo.addPlaneSurface([24], 24)    
        gmsh.model.geo.addPlaneSurface([25], 25)    
        gmsh.model.geo.addPlaneSurface([26], 26)
    end

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

    if second_beam
        gmsh.model.addPhysicalGroup(1, [21,22,23,24,25, 26, 27, 28, 29, 210, 211, 212], 12)
        gmsh.model.setPhysicalName(1, 12, "FreeEdges2")
        gmsh.model.addPhysicalGroup(2, [21, 22, 23, 25, 24], 23)
        gmsh.model.setPhysicalName(2, 23, "FreeAreas2")
        gmsh.model.addPhysicalGroup(2, [26], 24)
        gmsh.model.setPhysicalName(2, 24, "Dirichlet2")
        gmsh.model.addPhysicalGroup(3, [21], 32)
        gmsh.model.setPhysicalName(3, 32, "Beam2") 
    end

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
second_beam = true
x2 = 0.05
y2 = 1
z2 = 0.05
MeshGeneratorBeam(x1, y1, z1, x2, y2, z2, second_beam, 0.01)