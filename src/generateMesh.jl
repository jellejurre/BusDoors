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

    gmsh.model.geo.addPoint(w*0.25,   h*0.25,   0, meshSize, 5)
    gmsh.model.geo.addPoint(w*0.25,   h*0.75,   0, meshSize, 6)
    gmsh.model.geo.addPoint(w*0.75,   h*0.75,   0, meshSize, 7)
    gmsh.model.geo.addPoint(w*0.75,   h*0.25,   0, meshSize, 8)
    # Add lines
    gmsh.model.geo.addLine( 1,  2,  1)
    gmsh.model.geo.addLine( 2,  3,  2)
    gmsh.model.geo.addLine( 3,  4,  3)
    gmsh.model.geo.addLine( 1,  4,  4)

    gmsh.model.geo.addLine( 5,  6,  5)
    gmsh.model.geo.addLine( 6,  7,  6)
    gmsh.model.geo.addLine( 7,  8,  7)
    gmsh.model.geo.addLine( 5,  8,  8)
    # Construct curve loops and surfaces 
    gmsh.model.geo.addCurveLoop([1, 2, 3, -4], 1)
    gmsh.model.geo.addCurveLoop([5, 6, 7, -8], 2)

    gmsh.model.geo.addPlaneSurface([1, 2], 1)   #loop 1 with hole loop 2
    gmsh.model.geo.addPlaneSurface([2], 2)
    # Physical groups
    # def: addPhysicalGroup(dim, tags, tag = -1, name = "")
    #   Add a physical group of dimension `dim`, grouping the model entities with tags
    #   `tags`. Return the tag of the physical group, equal to `tag` if `tag` is
    #   positive, or a new tag if `tag` < 0. Set the name of the physical group if
    #   `name` is not empty.

    gmsh.model.addPhysicalGroup(1, [1,3], 2)
    gmsh.model.setPhysicalName(1, 2, "DirichletEdges")
    gmsh.model.addPhysicalGroup(1, [2,4, 5, 6, 7, 8], 3)
    gmsh.model.setPhysicalName(1, 3, "FreeEdges")
    gmsh.model.addPhysicalGroup(2, [1], 4)
    gmsh.model.setPhysicalName(2, 4, "FreeAreaOutside")
    gmsh.model.addPhysicalGroup(2, [2], 5)
    gmsh.model.setPhysicalName(2, 5, "FreeAreaInside")
    gmsh.model.geo.synchronize()
    # We can then generate a 2D mesh...
    gmsh.model.mesh.generate(2)
    # ... and save it to disk
    gmsh.write("geometry.msh")
    gmsh.finalize()
end

MeshGenerator(0.69, 1.93, 0.1)