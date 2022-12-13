using Gmsh
import Gmsh: gmsh
# 1 = aluminium frame, 2 = glass, 3 = rubber
# (x and y for rubber are thickness, not total width/height. No z value for rubber, goes from z1 to z2)
function MeshGenerator(x1, y1, z1, x2, y2, z2, z2offset, x3, y3, hinges_bool, hand_area_bool, meshSize)
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 1)
    gmsh.option.setNumber("Mesh.Algorithm", 6)
    gmsh.clear()
    gmsh.model.add("geometry")
    
    # low and high points of x2, y2 and z2
    x2l = 0.5*x1 - 0.5*x2
    x2h = 0.5*x1 + 0.5*x2
    y2l = 0.5*y1 - 0.5*y2
    y2h = 0.5*y1 + 0.5*y2
    z2l = 0.5*z1 - 0.5*z2 - z2offset
    z2h = 0.5*z1 + 0.5*z2 - z2offset

    x3l = 0.5*x1 - 0.5*x2 - x3
    x3h = 0.5*x1 + 0.5*x2 + x3
    y3l = 0.5*y1 - 0.5*y2 - y3
    y3h = 0.5*y1 + 0.5*y2 + y3
    z3l = 0
    z3h = z1

    # Tags: 1st digit = which physical group (aluminium = 1, glass = 2)
    # Add points
    gmsh.model.geo.addPoint(0,   0,   0, meshSize, 11)
    gmsh.model.geo.addPoint(0,   y1,   0, meshSize, 12)
    gmsh.model.geo.addPoint(x1,   y1,   0, meshSize, 13)
    gmsh.model.geo.addPoint(x1,   0,   0, meshSize, 14)
    gmsh.model.geo.addPoint(0,   0,   z1, meshSize, 15)
    gmsh.model.geo.addPoint(0,   y1,   z1, meshSize, 16)
    gmsh.model.geo.addPoint(x1,   y1,   z1, meshSize, 17)
    gmsh.model.geo.addPoint(x1,   0,   z1, meshSize, 18)

    gmsh.model.geo.addPoint(x2l,   y2l,   z2l, meshSize, 21)
    gmsh.model.geo.addPoint(x2l,   y2h,   z2l, meshSize, 22)
    gmsh.model.geo.addPoint(x2h,   y2h,   z2l, meshSize, 23)
    gmsh.model.geo.addPoint(x2h,   y2l,   z2l, meshSize, 24)
    gmsh.model.geo.addPoint(x2l,   y2l,   z2h, meshSize, 25)
    gmsh.model.geo.addPoint(x2l,   y2h,   z2h, meshSize, 26)
    gmsh.model.geo.addPoint(x2h,   y2h,   z2h, meshSize, 27)
    gmsh.model.geo.addPoint(x2h,   y2l,   z2h, meshSize, 28)

    gmsh.model.geo.addPoint(x3l,   y3l,   z3l, meshSize, 31)
    gmsh.model.geo.addPoint(x3l,   y3h,   z3l, meshSize, 32)
    gmsh.model.geo.addPoint(x3h,   y3h,   z3l, meshSize, 33)
    gmsh.model.geo.addPoint(x3h,   y3l,   z3l, meshSize, 34)
    gmsh.model.geo.addPoint(x3l,   y3l,   z3h, meshSize, 35)
    gmsh.model.geo.addPoint(x3l,   y3h,   z3h, meshSize, 36)
    gmsh.model.geo.addPoint(x3h,   y3h,   z3h, meshSize, 37)
    gmsh.model.geo.addPoint(x3h,   y3l,   z3h, meshSize, 38)
    
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

    gmsh.model.geo.addLine( 31,  32,  31)
    gmsh.model.geo.addLine( 32,  33,  32)
    gmsh.model.geo.addLine( 33,  34,  33)
    gmsh.model.geo.addLine( 31,  34,  34)
    gmsh.model.geo.addLine( 35,  36,  35)
    gmsh.model.geo.addLine( 36,  37,  36)
    gmsh.model.geo.addLine( 37,  38,  37)
    gmsh.model.geo.addLine( 35,  38,  38)
    gmsh.model.geo.addLine( 35,  31,  39)
    gmsh.model.geo.addLine( 36,  32,  310)
    gmsh.model.geo.addLine( 37,  33,  311)
    gmsh.model.geo.addLine( 38,  34,  312)

    # funnel lines
    gmsh.model.geo.addLine( 31,  21,  41)
    gmsh.model.geo.addLine( 32,  22,  42)
    gmsh.model.geo.addLine( 33,  23,  43)
    gmsh.model.geo.addLine( 34,  24,  44)
    gmsh.model.geo.addLine( 35,  25,  45)
    gmsh.model.geo.addLine( 36,  26,  46)
    gmsh.model.geo.addLine( 37,  27,  47)
    gmsh.model.geo.addLine( 38,  28,  48)

    # Construct curve loops and surfaces    (watch out for line direction)
    gmsh.model.geo.addCurveLoop([11, 12, 13, -14], 11)
    gmsh.model.geo.addCurveLoop([15, 16, 17, -18], 12)
    gmsh.model.geo.addCurveLoop([11, -110, -15, 19], 13)
    gmsh.model.geo.addCurveLoop([-13, -111, 17, 112], 14)
    gmsh.model.geo.addCurveLoop([-14, -19, 18, 112], 15)
    gmsh.model.geo.addCurveLoop([111, -12, -110, 16], 16)
    gmsh.model.geo.addCurveLoop([21, 22, 23, -24], 21)
    gmsh.model.geo.addCurveLoop([25, 26, 27, -28], 22)
    gmsh.model.geo.addCurveLoop([21, -210, -25, 29], 23)
    gmsh.model.geo.addCurveLoop([-23, -211, 27, 212], 24)
    gmsh.model.geo.addCurveLoop([-24, -29, 28, 212], 25)
    gmsh.model.geo.addCurveLoop([211, -22, -210, 26], 26)

    gmsh.model.geo.addCurveLoop([31, -310, -35, 39], 31)
    gmsh.model.geo.addCurveLoop([311, -32, -310, 36], 32)
    gmsh.model.geo.addCurveLoop([-33, -311, 37, 312], 33)
    gmsh.model.geo.addCurveLoop([-34, -39, 38, 312], 34)
    gmsh.model.geo.addCurveLoop([46, -25, -45, 35], 35)
    gmsh.model.geo.addCurveLoop([47, -26, -46, 36], 36)
    gmsh.model.geo.addCurveLoop([48, -27, -47, 37], 37)
    gmsh.model.geo.addCurveLoop([45, 28, -48, -38], 38)
    gmsh.model.geo.addCurveLoop([42, -21, -41, 31], 39)
    gmsh.model.geo.addCurveLoop([43, -22, -42, 32], 310)
    gmsh.model.geo.addCurveLoop([44, -23, -43, 33], 311)
    gmsh.model.geo.addCurveLoop([41, 24, -44, -34], 312)
    gmsh.model.geo.addCurveLoop([46, 210, -42, -310], 313)
    gmsh.model.geo.addCurveLoop([47, 211, -43, -311], 314)
    gmsh.model.geo.addCurveLoop([48, 212, -44, -312], 315)
    gmsh.model.geo.addCurveLoop([45, 29, -41, -39], 316)

    # loops for the holes in aluminium frame
    gmsh.model.geo.addCurveLoop([31, 32, 33, -34], 317)
    gmsh.model.geo.addCurveLoop([35, 36, 37, -38], 318)



    gmsh.model.geo.addPlaneSurface([11, 317], 11)   #loop with hole loop inside
    gmsh.model.geo.addPlaneSurface([12, 318], 12)
    gmsh.model.geo.addPlaneSurface([13], 13)
    gmsh.model.geo.addPlaneSurface([14], 14)
    gmsh.model.geo.addPlaneSurface([15], 15)
    gmsh.model.geo.addPlaneSurface([16], 16)
    gmsh.model.geo.addPlaneSurface([21], 21)
    gmsh.model.geo.addPlaneSurface([22], 22)    
    gmsh.model.geo.addPlaneSurface([23], 23)    
    gmsh.model.geo.addPlaneSurface([24], 24)    
    gmsh.model.geo.addPlaneSurface([25], 25)    
    gmsh.model.geo.addPlaneSurface([26], 26)
    gmsh.model.geo.addPlaneSurface([31], 31)
    gmsh.model.geo.addPlaneSurface([32], 32)    
    gmsh.model.geo.addPlaneSurface([33], 33)    
    gmsh.model.geo.addPlaneSurface([34], 34)    
    gmsh.model.geo.addPlaneSurface([35], 35)    
    gmsh.model.geo.addPlaneSurface([36], 36)
    gmsh.model.geo.addPlaneSurface([37], 37)
    gmsh.model.geo.addPlaneSurface([38], 38)
    gmsh.model.geo.addPlaneSurface([39], 39)
    gmsh.model.geo.addPlaneSurface([310], 310)
    gmsh.model.geo.addPlaneSurface([311], 311)
    gmsh.model.geo.addPlaneSurface([312], 312)
    gmsh.model.geo.addPlaneSurface([313], 313)
    gmsh.model.geo.addPlaneSurface([314], 314)
    gmsh.model.geo.addPlaneSurface([315], 315)
    gmsh.model.geo.addPlaneSurface([316], 316)

    # outer volume aluminium
    gmsh.model.geo.addSurfaceLoop([11, 12, 13, 14, 15, 16], 11)
    gmsh.model.geo.addSurfaceLoop([31, 34, 33, 32], 12)     #only remove the sides, already a hole in surface 11, 12
    gmsh.model.geo.addVolume([11, 12], 11)
    
    # inner volume glass
    gmsh.model.geo.addSurfaceLoop([21, 22, 23, 24, 25, 26], 21)
    gmsh.model.geo.addVolume([21], 21)

    # volume rubber
    gmsh.model.geo.addSurfaceLoop([31, 35, 39, 23, 313, 316], 31)
    gmsh.model.geo.addVolume([31], 31)
    gmsh.model.geo.addSurfaceLoop([32, 36, 310, 26, 313, 314], 32)
    gmsh.model.geo.addVolume([32], 32)
    gmsh.model.geo.addSurfaceLoop([33, 37, 311, 24, 314, 315], 33)
    gmsh.model.geo.addVolume([33], 33)
    gmsh.model.geo.addSurfaceLoop([34, 38, 312, 25, 315, 316], 34)
    gmsh.model.geo.addVolume([34], 34)

    # cylinder (top left)
    if hinges_bool
        # location center of cylinder
        xcyl1 = x1-0.04107
        ycyl1 = y1-0.091677
        zcyl1 = 0
        cylheight1 = -0.05  # -, since it pokes out to negative z
        r = 0.015458
        gmsh.model.geo.addPoint(xcyl1,   ycyl1,   zcyl1, meshSize, 1000)
        gmsh.model.geo.addPoint(xcyl1 + r,   ycyl1,   zcyl1, meshSize, 1001)
        gmsh.model.geo.addPoint(xcyl1,   ycyl1 - r,   zcyl1, meshSize, 1002)
        gmsh.model.geo.addPoint(xcyl1 - r,   ycyl1,   zcyl1, meshSize, 1003)
        gmsh.model.geo.addPoint(xcyl1,   ycyl1 + r,   zcyl1, meshSize, 1004)
        # gmsh.model.geo.addCircleArc(startTag, centerTag, endTag, tag = -1, nx = 0., ny = 0., nz = 0.), must have arc < pi
        gmsh.model.geo.addCircleArc(1001, 1000, 1002, 1000)
        gmsh.model.geo.addCircleArc(1002, 1000, 1003, 1001)
        gmsh.model.geo.addCircleArc(1003, 1000, 1004, 1002)
        gmsh.model.geo.addCircleArc(1004, 1000, 1001, 1003)
        gmsh.model.geo.addCurveLoop([1000, 1001, 1002, 1003], 1000)
        gmsh.model.geo.addPlaneSurface([1000], 1000)
        # extrude returns array of pairs (dim, tag) for volumes
        cyltags = gmsh.model.geo.extrude([(2, 1000)], 0, 0, cylheight1)
        # filter away 2d planes, keeping 3d volumes
        cyl_ceilingtag = cyltags[1][2]
        cyl_volume_tag = filter(tuple -> tuple[1] == 3, cyltags)
        cyl_volume_tag = map(tuple -> tuple[2], cyl_volume_tag)[1]
        free_cyl_tags = symdiff(map(tuple -> tuple[2], cyltags), [cyl_ceilingtag, cyl_volume_tag])

        
        # hinge box (bottom right)
        # dimensions box:
        xbox1 = 0.05
        ybox1 = 0.05

        # locations points:
        xbox1l = 0.05
        ybox1l = 0.05
        zbox1 = 0
        boxheight1 = -0.05 # -, since it pokes out to negative z
        gmsh.model.geo.addPoint(xbox1l,   ybox1l,   zbox1, meshSize, 1100)
        gmsh.model.geo.addPoint(xbox1l + xbox1,   ybox1l,   zbox1, meshSize, 1101)
        gmsh.model.geo.addPoint(xbox1l,   ybox1l + ybox1,   zbox1, meshSize, 1102)
        gmsh.model.geo.addPoint(xbox1l + xbox1,   ybox1l + ybox1,   zbox1, meshSize, 1103)

        gmsh.model.geo.addLine( 1100,  1101,  1100)
        gmsh.model.geo.addLine( 1101,  1103,  1101)
        gmsh.model.geo.addLine( 1103,  1102,  1102)
        gmsh.model.geo.addLine( 1102,  1100,  1103)
        gmsh.model.geo.addCurveLoop([1100, 1101, 1102, 1103], 1100)
        gmsh.model.geo.addPlaneSurface([1100], 1100)
        # extrude returns array of pairs (dim, tag) for volumes
        boxtags = gmsh.model.geo.extrude([(2, 1100)], 0, 0, boxheight1)
        # filter away 2d planes, keeping 3d volumes
        box_ceilingtag = boxtags[1][2]
        box_volume_tag = filter(tuple -> tuple[1] == 3, boxtags)
        box_volume_tag = map(tuple -> tuple[2], box_volume_tag)[1]
        free_box_tags = symdiff(map(tuple -> tuple[2], boxtags), [box_ceilingtag, box_volume_tag])

        gmsh.model.addPhysicalGroup(3, [cyl_volume_tag, box_volume_tag], 1000)
        gmsh.model.setPhysicalName(3, 1000, "Hinges")
        gmsh.model.addPhysicalGroup(2, [cyl_ceilingtag, box_ceilingtag], 1001)
        gmsh.model.setPhysicalName(2, 1001, "Hinge ceilings")
        gmsh.model.addPhysicalGroup(2, append!(free_box_tags, free_cyl_tags), 1002)
        gmsh.model.setPhysicalName(2, 1002, "Hinge sides")
    end

    if hand_area_bool
        # location center of hand
        xhand = x1*0.5
        yhand = y1*0.5
        zhand = 0
        r = 0.0315
        gmsh.model.geo.addPoint(xhand,   yhand,   zhand, meshSize, 1200)
        gmsh.model.geo.addPoint(xhand + r,   yhand,   zhand, meshSize, 1201)
        gmsh.model.geo.addPoint(xhand,   yhand - r,   zhand, meshSize, 1202)
        gmsh.model.geo.addPoint(xhand - r,   yhand,   zhand, meshSize, 1203)
        gmsh.model.geo.addPoint(xhand,   yhand + r,   zhand, meshSize, 1204)
        
        # gmsh.model.geo.addCircleArc(startTag, centerTag, endTag, tag = -1, nx = 0., ny = 0., nz = 0.), must have arc < pi
        gmsh.model.geo.addCircleArc(1201, 1200, 1202, 1200)
        gmsh.model.geo.addCircleArc(1202, 1200, 1203, 1201)
        gmsh.model.geo.addCircleArc(1203, 1200, 1204, 1202)
        gmsh.model.geo.addCircleArc(1204, 1200, 1201, 1203)
        gmsh.model.geo.addCurveLoop([1200, 1201, 1202, 1203], 1200)
        gmsh.model.geo.addPlaneSurface([1200], 1200)

        gmsh.model.addPhysicalGroup(2, [1200], 1200)
        gmsh.model.setPhysicalName(2, 1200, "Hand")
    end

    # Physical groups
    # def: addPhysicalGroup(dim, tags, tag = -1, name = "")
    #   Add a physical group of dimension `dim`, grouping the model entities with tags
    #   `tags`. Return the tag of the physical group, equal to `tag` if `tag` is
    #   positive, or a new tag if `tag` < 0. Set the name of the physical group if
    #   `name` is not empty.

    gmsh.model.addPhysicalGroup(1, [11,12,13,14,15, 16, 17, 18, 19, 110, 111, 112, 21, 22, 23, 24, 25, 26, 27, 28, 29, 210, 211, 212, 31, 32, 33, 34, 35, 36, 37, 38, 39, 310, 311, 312, 1000, 1001, 1002, 1003], 11)
    gmsh.model.setPhysicalName(1, 11, "FreeEdges")
    gmsh.model.addPhysicalGroup(2, [11, 12, 13, 14, 15, 16, 21, 22, 23, 24, 25, 26, 31, 32, 33, 34, 35, 36, 37, 38, 39, 310, 311, 312, 313, 314, 315, 316], 2)
    gmsh.model.setPhysicalName(2, 2, "FreeAreas")
    gmsh.model.addPhysicalGroup(3, [11], 31)
    gmsh.model.setPhysicalName(3, 31, "Aluminium")
    gmsh.model.addPhysicalGroup(3, [21], 32)
    gmsh.model.setPhysicalName(3, 32, "Glass")
    gmsh.model.addPhysicalGroup(3, [31, 32, 33, 34], 33)
    gmsh.model.setPhysicalName(3, 33, "Rubber")

    gmsh.model.geo.synchronize()
    # We can then generate a 3D mesh...
    gmsh.model.mesh.generate(3)
    # ... and save it to disk
    gmsh.write("geometry.msh")
    gmsh.finalize()
end

x1 = 0.6610
y1 = 2.02366
z1 = 0.0323
x2 = 0.4855
y2 = 1.764
z2 = 0.0023
z2offset = z1*0.5
x3 = 0.001531
y3 = 0.005081

# sizes and location of the following parts should be changed in the function MeshGenerator above
# booleans for enabling/disabling certain parts of the mesh:
# Physical group: Hinges, Hinge ceilings, Hinge sides
hinges_bool = true
# Physical group: Hand
hand_area_bool = true

MeshGenerator(x1, y1, z1, x2, y2, z2, z2offset, x3, y3, hinges_bool, hand_area_bool, 0.01)