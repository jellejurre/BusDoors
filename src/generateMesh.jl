using Gmsh
import Gmsh: gmsh

# At the bottom of the code the input values can be changed for generating the model

# Function to generate the parametric door and optionally two hinges and a circular area for a pushing hand.
# aluminium: tuple of x, y, z dimensions of the outer aluminium frame
# glass: tuple of x, y, z dimensions of the inner glass plane
# rubber: tuple of x, y -> the width of the rubber between the aluminium and glass on the x axis and y axis
# hinges_bool: if true, add two hinges (a cylinder and a box) to the model
# hand_area_bool: if true, add a circular area for a pushing hand.
# cyl_measurements: tuple of x, y, z, cylheight, rcyl -> the location of the center point of the cylinder, the height of the cylinder and the radius of the circle
# box_measurements: tuple of xbox, ybox, boxheight, x, y, z -> the dimensions of the box in x, y and its height. Then the location of the point of the corner closest to (0,0,0)
# hand_measurements: tuple of x, y, z, r -> location of the hand and its radius
# meshSize: target mesh sizes

function MeshGenerator(aluminium, glass, rubber, hinges_bool, hand_area_bool, cyl_measurements, box_measurements, hand_measurements, meshSize)
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 1)
    gmsh.option.setNumber("Mesh.Algorithm", 6)
    gmsh.clear()
    gmsh.model.add("geometry")

    aluminium_x = aluminium[1]
    aluminium_y = aluminium[2]
    aluminium_z = aluminium[3]

    glass_x = glass[1]
    glass_y = glass[2]
    glass_z = glass[3]

    # the width of the rubber between the aluminium and glass on the x axis and y axis
    rubber_x = rubber[1]
    rubber_y = rubber[2]
    
    # low and high points of glass_x, glass_y and glass_z
    glass_xl = 0.5*aluminium_x - 0.5*glass_x
    glass_xh = 0.5*aluminium_x + 0.5*glass_x
    glass_yl = 0.5*aluminium_y - 0.5*glass_y
    glass_yh = 0.5*aluminium_y + 0.5*glass_y
    glass_zl = aluminium_z - glass_z
    glass_zh = aluminium_z

    rubber_xl = 0.5*aluminium_x - 0.5*glass_x - rubber_x
    rubber_xh = 0.5*aluminium_x + 0.5*glass_x + rubber_x
    rubber_yl = 0.5*aluminium_y - 0.5*glass_y - rubber_y
    rubber_yh = 0.5*aluminium_y + 0.5*glass_y + rubber_y
    z3l = 0
    z3h = aluminium_z

    # Tags: 1st digit = which physical group (aluminium = 1, glass = 2)
    # Add points
    gmsh.model.geo.addPoint(0,   0,   0, meshSize, 11)
    gmsh.model.geo.addPoint(0,   aluminium_y,   0, meshSize, 12)
    gmsh.model.geo.addPoint(aluminium_x,   aluminium_y,   0, meshSize, 13)
    gmsh.model.geo.addPoint(aluminium_x,   0,   0, meshSize, 14)
    gmsh.model.geo.addPoint(0,   0,   aluminium_z, meshSize, 15)
    gmsh.model.geo.addPoint(0,   aluminium_y,   aluminium_z, meshSize, 16)
    gmsh.model.geo.addPoint(aluminium_x,   aluminium_y,   aluminium_z, meshSize, 17)
    gmsh.model.geo.addPoint(aluminium_x,   0,   aluminium_z, meshSize, 18)

    gmsh.model.geo.addPoint(glass_xl,   glass_yl,   glass_zl, meshSize, 21)
    gmsh.model.geo.addPoint(glass_xl,   glass_yh,   glass_zl, meshSize, 22)
    gmsh.model.geo.addPoint(glass_xh,   glass_yh,   glass_zl, meshSize, 23)
    gmsh.model.geo.addPoint(glass_xh,   glass_yl,   glass_zl, meshSize, 24)
    gmsh.model.geo.addPoint(glass_xl,   glass_yl,   glass_zh, meshSize, 25)
    gmsh.model.geo.addPoint(glass_xl,   glass_yh,   glass_zh, meshSize, 26)
    gmsh.model.geo.addPoint(glass_xh,   glass_yh,   glass_zh, meshSize, 27)
    gmsh.model.geo.addPoint(glass_xh,   glass_yl,   glass_zh, meshSize, 28)

    gmsh.model.geo.addPoint(rubber_xl,   rubber_yl,   z3l, meshSize, 31)
    gmsh.model.geo.addPoint(rubber_xl,   rubber_yh,   z3l, meshSize, 32)
    gmsh.model.geo.addPoint(rubber_xh,   rubber_yh,   z3l, meshSize, 33)
    gmsh.model.geo.addPoint(rubber_xh,   rubber_yl,   z3l, meshSize, 34)
    gmsh.model.geo.addPoint(rubber_xl,   rubber_yl,   z3h, meshSize, 35)
    gmsh.model.geo.addPoint(rubber_xl,   rubber_yh,   z3h, meshSize, 36)
    gmsh.model.geo.addPoint(rubber_xh,   rubber_yh,   z3h, meshSize, 37)
    gmsh.model.geo.addPoint(rubber_xh,   rubber_yl,   z3h, meshSize, 38)
    
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
    if !hinges_bool
        gmsh.model.geo.addSurfaceLoop([11, 12, 13, 14, 15, 16], 11)
        gmsh.model.geo.addSurfaceLoop([31, 34, 33, 32], 12)     #only remove the sides, already a hole in surface 11, 12
        gmsh.model.geo.addVolume([11, 12], 11)
    end
    if hand_area_bool
        # location center of hand
        xhand = hand_measurements[1]
        yhand = hand_measurements[2]
        zhand = hand_measurements[3]
        rhand = hand_measurements[4]
        gmsh.model.geo.addPoint(xhand,   yhand,   zhand, meshSize, 1200)
        gmsh.model.geo.addPoint(xhand + rhand,   yhand,   zhand, meshSize, 1201)
        gmsh.model.geo.addPoint(xhand,   yhand - rhand,   zhand, meshSize, 1202)
        gmsh.model.geo.addPoint(xhand - rhand,   yhand,   zhand, meshSize, 1203)
        gmsh.model.geo.addPoint(xhand,   yhand + rhand,   zhand, meshSize, 1204)
        
        # gmsh.model.geo.addCircleArc(startTag, centerTag, endTag, tag = -1, nx = 0., ny = 0., nz = 0.), must have arc < pi
        gmsh.model.geo.addCircleArc(1201, 1200, 1202, 1200)
        gmsh.model.geo.addCircleArc(1202, 1200, 1203, 1201)
        gmsh.model.geo.addCircleArc(1203, 1200, 1204, 1202)
        gmsh.model.geo.addCircleArc(1204, 1200, 1201, 1203)
        gmsh.model.geo.addCurveLoop([1200, 1201, 1202, 1203], 1200)
        gmsh.model.geo.addPlaneSurface([1200], 1200)

        gmsh.model.addPhysicalGroup(2, [1200], 1200)
        gmsh.model.setPhysicalName(2, 1200, "Hand")
        gmsh.model.addPhysicalGroup(1, [1200, 1201, 1202, 1203], 1201)
        gmsh.model.setPhysicalName(1, 1201, "Hand edges")
    end
    if hand_area_bool
        gmsh.model.geo.addPlaneSurface([21, 1200], 1202) # glass area - hand

        # inner volume glass
        gmsh.model.geo.addSurfaceLoop([1202, 1200, 22, 23, 24, 25, 26], 21)
        gmsh.model.geo.addVolume([21], 21)
    else
    # inner volume glass
    gmsh.model.geo.addSurfaceLoop([21, 22, 23, 24, 25, 26], 21)
    gmsh.model.geo.addVolume([21], 21)
    end

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
        xcyl = cyl_measurements[1]
        ycyl = cyl_measurements[2]
        zcyl = cyl_measurements[3]
        cylheight = cyl_measurements[4]
        rcyl = cyl_measurements[5]
        gmsh.model.geo.addPoint(xcyl,   ycyl,   zcyl, meshSize, 1000)
        gmsh.model.geo.addPoint(xcyl + rcyl,   ycyl,   zcyl, meshSize, 1001)
        gmsh.model.geo.addPoint(xcyl,   ycyl - rcyl,   zcyl, meshSize, 1002)
        gmsh.model.geo.addPoint(xcyl - rcyl,   ycyl,   zcyl, meshSize, 1003)
        gmsh.model.geo.addPoint(xcyl,   ycyl + rcyl,   zcyl, meshSize, 1004)
        # gmsh.model.geo.addCircleArc(startTag, centerTag, endTag, tag = -1, nx = 0., ny = 0., nz = 0.), must have arc < pi
        gmsh.model.geo.addCircleArc(1001, 1000, 1002, 1000)
        gmsh.model.geo.addCircleArc(1002, 1000, 1003, 1001)
        gmsh.model.geo.addCircleArc(1003, 1000, 1004, 1002)
        gmsh.model.geo.addCircleArc(1004, 1000, 1001, 1003)
        gmsh.model.geo.addCurveLoop([1000, 1001, 1002, 1003], 1000)
        gmsh.model.geo.addPlaneSurface([1000], 1000)
        # extrude returns array of pairs (dim, tag) for volumes
        cyltags = gmsh.model.geo.extrude([(2, 1000)], 0, 0, cylheight)
        # filter away 2d planes, keeping 3d volumes
        cyl_ceilingtag = cyltags[1][2]
        cyl_volume_tag = filter(tuple -> tuple[1] == 3, cyltags)
        cyl_volume_tag = map(tuple -> tuple[2], cyl_volume_tag)[1]
        free_cyl_tags = symdiff(map(tuple -> tuple[2], cyltags), [cyl_ceilingtag, cyl_volume_tag])

        
        # hinge box (bottom right)
        # dimensions box:
        xbox = box_measurements[1]
        ybox = box_measurements[2]
        boxheight = box_measurements[3]
        # locations points:
        xboxl = box_measurements[4]
        yboxl = box_measurements[5]
        zbox = box_measurements[6]

        gmsh.model.geo.addPoint(xboxl,   yboxl,   zbox, meshSize, 1100)
        gmsh.model.geo.addPoint(xboxl + xbox,   yboxl,   zbox, meshSize, 1101)
        gmsh.model.geo.addPoint(xboxl,   yboxl + ybox,   zbox, meshSize, 1102)
        gmsh.model.geo.addPoint(xboxl + xbox,   yboxl + ybox,   zbox, meshSize, 1103)

        gmsh.model.geo.addLine( 1100,  1101,  1100)
        gmsh.model.geo.addLine( 1101,  1103,  1101)
        gmsh.model.geo.addLine( 1103,  1102,  1102)
        gmsh.model.geo.addLine( 1102,  1100,  1103)
        gmsh.model.geo.addCurveLoop([1100, 1101, 1102, 1103], 1100)
        gmsh.model.geo.addPlaneSurface([1100], 1100)
        # extrude returns array of pairs (dim, tag) for volumes
        boxtags = gmsh.model.geo.extrude([(2, 1100)], 0, 0, boxheight)
        # filter away 2d planes, keeping 3d volumes
        box_ceilingtag = boxtags[1][2]
        box_volume_tag = filter(tuple -> tuple[1] == 3, boxtags)
        box_volume_tag = map(tuple -> tuple[2], box_volume_tag)[1]
        free_box_tags = symdiff(map(tuple -> tuple[2], boxtags), [box_ceilingtag, box_volume_tag])

        gmsh.model.geo.addPlaneSurface([11, 317, 1000, 1100], 1101)   #loop with hole loop inside and holes for hinges
        gmsh.model.geo.addSurfaceLoop([1101, 12, 13, 14, 15, 16, 1000, 1100], 11)
        gmsh.model.geo.addSurfaceLoop([31, 34, 33, 32], 12)     #only remove the sides, already a hole in surface 11, 12
        gmsh.model.geo.addVolume([11, 12], 11)

        gmsh.model.addPhysicalGroup(1, [1000, 1001, 1002, 1003, 1100, 1101, 1102, 1103], 1003)
        gmsh.model.setPhysicalName(1, 1003, "Hinge edges")

        gmsh.model.addPhysicalGroup(3, [cyl_volume_tag], 1000)
        gmsh.model.setPhysicalName(3, 1000, "Hinge top")
        gmsh.model.addPhysicalGroup(2, [cyl_ceilingtag], 1001)
        gmsh.model.setPhysicalName(2, 1001, "Hinge ceiling top")
        gmsh.model.addPhysicalGroup(2, free_cyl_tags, 1002)
        gmsh.model.setPhysicalName(2, 1002, "Hinge sides top")

        gmsh.model.addPhysicalGroup(3, [box_volume_tag], 1100)
        gmsh.model.setPhysicalName(3, 1100, "Hinge bottom")
        gmsh.model.addPhysicalGroup(2, [box_ceilingtag], 1101)
        gmsh.model.setPhysicalName(2, 1101, "Hinge ceiling bottom")
        gmsh.model.addPhysicalGroup(2, free_box_tags, 1102)
        gmsh.model.setPhysicalName(2, 1102, "Hinge sides bottom")
    end



    # Physical groups
    # def: addPhysicalGroup(dim, tags, tag = -1, name = "")
    #   Add a physical group of dimension `dim`, grouping the model entities with tags
    #   `tags`. Return the tag of the physical group, equal to `tag` if `tag` is
    #   positive, or a new tag if `tag` < 0. Set the name of the physical group if
    #   `name` is not empty.

    gmsh.model.addPhysicalGroup(1, [11,12,13,14,15, 16, 17, 18, 19, 110, 111, 112, 21, 22, 23, 24, 25, 26, 27, 28, 29, 210, 211, 212, 31, 32, 33, 34, 35, 36, 37, 38, 39, 310, 311, 312], 11)
    gmsh.model.setPhysicalName(1, 11, "FreeEdges")
    if hand_area_bool
        if hinges_bool
            gmsh.model.addPhysicalGroup(2, [1101, 12, 1202, 22, 23, 24, 25, 26, 31, 32, 33, 34, 35, 36, 37, 38, 39, 310, 311, 312, 313, 314, 315, 316], 2)
            gmsh.model.setPhysicalName(2, 2, "FreeAreas")
        else
            gmsh.model.addPhysicalGroup(2, [11, 12, 1202, 22, 23, 24, 25, 26, 31, 32, 33, 34, 35, 36, 37, 38, 39, 310, 311, 312, 313, 314, 315, 316], 2)
            gmsh.model.setPhysicalName(2, 2, "FreeAreas")
        end
    else
        if hinges_bool
            gmsh.model.addPhysicalGroup(2, [1101, 12, 21, 22, 23, 24, 25, 26, 31, 32, 33, 34, 35, 36, 37, 38, 39, 310, 311, 312, 313, 314, 315, 316], 2)
            gmsh.model.setPhysicalName(2, 2, "FreeAreas")
        else
            gmsh.model.addPhysicalGroup(2, [11, 12, 21, 22, 23, 24, 25, 26, 31, 32, 33, 34, 35, 36, 37, 38, 39, 310, 311, 312, 313, 314, 315, 316], 2)
            gmsh.model.setPhysicalName(2, 2, "FreeAreas")
        end
    end
    gmsh.model.addPhysicalGroup(2, [16], 21)
    gmsh.model.setPhysicalName(2, 21, "Area top")
    gmsh.model.addPhysicalGroup(2, [15], 22)
    gmsh.model.setPhysicalName(2, 22, "Area bottom")
    gmsh.model.addPhysicalGroup(2, [13], 23)
    gmsh.model.setPhysicalName(2, 23, "Area left")
    gmsh.model.addPhysicalGroup(2, [14], 24)
    gmsh.model.setPhysicalName(2, 24, "Area right")
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

# dimensions of the door
aluminium_x = 0.6610
aluminium_y = 2.02366
aluminium_z = 0.0323
aluminium = (aluminium_x, aluminium_y, aluminium_z)

glass_x = 0.4855
glass_y = 1.764
glass_z = 0.0023
glass = (glass_x, glass_y, glass_z)

# the width of the rubber between the aluminium and glass on the x axis and y axis
rubber_x = 0.001531
rubber_y = 0.005081
rubber = (rubber_x, rubber_y)

# hinge cylinder (top left)
# location center of cylinder
xcyl = aluminium_x-0.04107
ycyl = aluminium_y-0.091677
zcyl = 0
# dimensions cylinder:
cylheight = -0.05  # -, since it pokes out to negative z
rcyl = 0.015458
cyl_measurements = (xcyl, ycyl, zcyl, cylheight, rcyl)

# hinge box (bottom right)
# dimensions box:
xbox = 0.05
ybox = 0.05
boxheight = -0.05 # -, since it pokes out to negative z
# location of point of corner closest to (0,0,0):
xboxl = 0.05
yboxl = 0.05
zbox = 0
box_measurements = (xbox, ybox, boxheight, xboxl, yboxl, zbox)

# location center of hand
xhand = aluminium_x*0.5
yhand = aluminium_y*0.5
zhand = aluminium_z - glass_z
rhand = 0.1  # radius of the circle representing the hand
hand_measurements = (xhand, yhand, zhand, rhand)

# sizes and location of the following parts should be changed in the function MeshGenerator above
# booleans for enabling/disabling certain parts of the mesh:
# Physical group: Hinge top/bottom, Hinge ceiling top/bottom, Hinge sides top/bottom, Hinge edges
hinges_bool = true
# Physical group: Hand, Hand edges
hand_area_bool = true

MeshGenerator(aluminium, glass, rubber, hinges_bool, hand_area_bool, cyl_measurements, box_measurements, hand_measurements, 0.01)