using Gmsh
import Gmsh: gmsh
function MeshGenerator(w, h, meshSize)
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 1)
    gmsh.option.setNumber("Mesh.Algorithm", 6)
    gmsh.clear()
    gmsh.model.add("geometry")
    # Add points
    gmsh.model.geo.addPoint(0,   0,   0, meshSize, 1)
    gmsh.model.geo.addPoint(0,   h,   0, meshSize, 2)
    gmsh.model.geo.addPoint(w,   h,   0, meshSize, 3)
    gmsh.model.geo.addPoint(w,   0,   0, meshSize, 4)
    # Add lines
    gmsh.model.geo.addLine( 1,  2,  1)
    gmsh.model.geo.addLine( 2,  3,  2)
    gmsh.model.geo.addLine( 3,  4,  3)
    gmsh.model.geo.addLine( 1,  4,  4)
    # Construct curve loops and surfaces 
    gmsh.model.geo.addCurveLoop([1, 2, 3, -4], 1)
    gmsh.model.geo.addPlaneSurface([1], 1)
    # Physical groups
    gmsh.model.addPhysicalGroup(1, [1,3], 2)
    gmsh.model.setPhysicalName(1, 2, "DirichletEdges")
    gmsh.model.addPhysicalGroup(1, [2,4], 3)
    gmsh.model.setPhysicalName(1, 3, "FreeEdges")
    gmsh.model.addPhysicalGroup(2, [1], 4)
    gmsh.model.setPhysicalName(2, 4, "FreeArea")
    gmsh.model.geo.synchronize()
    # We can then generate a 2D mesh...
    gmsh.model.mesh.generate(2)
    # ... and save it to disk
    gmsh.write("geometry.msh")
    gmsh.finalize()
end

MeshGenerator(0.69, 1.93, 0.01)