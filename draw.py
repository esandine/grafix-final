from display import *
from matrix import *
from math import *
from gmath import *
import random

#coeffs in form of r_amb, r_diff, r_spec, g.., b..
def gen_color(matrix, point, lights, normal):
    r = 0
    g = 0
    b = 0
    for light in lights:
        Iamb = gen_iamb(matrix, point, light)
        Idiff = gen_idiff(matrix, point, light, normal)
        Ispec = gen_ispec(matrix, point, light, normal)
        coeffs = light[1]['constants']
        r += coeffs['red'][0]*Iamb[0]+coeffs['red'][1]*Idiff[0]+coeffs['red'][2]*Ispec[0]
        g += coeffs['green'][0]*Iamb[1]+coeffs['green'][1]*Idiff[1]+coeffs['green'][2]*Ispec[1]
        b += coeffs['blue'][0]*Iamb[2]+coeffs['blue'][1]*Idiff[2]+coeffs['blue'][2]*Ispec[2]
    if r > 255:
        r=255
    if g > 255:
        g=255
    if b > 255:
        b=255
    if r < 0:
        r=0
    if g < 0:
        g=0
    if b < 0:
        b=0
    
    return r,g,b

def gen_iamb(matrix, point, light):
    return light[1]['color']

def normalize(xyz):
    m = mag(xyz)
    if m == 0:
        return 0, 0, 0
    else:
        return mult(xyz, 1/m)

def mag(xyz):
    x = xyz[0]
    y = xyz[1]
    z = xyz[2]
    return math.sqrt(x*x+y*y+z*z)

def add(v0, v1):
    return v0[0]+v1[0], v0[1]+v1[1], v0[2]+v1[2]

def mult(v0, s):
    return v0[0]*s, v0[1]*s, v0[2]*s

def sub(v0, v1):
    return add(v0, mult(v1,-1))

def dot(v0, v1):
    return v0[0]*v1[0]+v0[1]*v1[1]+v0[2]*v1[2]

def costheta(v0, v1):#does dot product stuff
    v0 = normalize(v0)
    v1 = normalize(v1)
    return v0[0]*v1[0]+v0[1]*v1[1]+v0[2]*v1[2]

def gen_idiff(matrix, point, light, normal):
    ret = [0,0,0]
    l = light[1]['location']
    ctheta = costheta(l, normal)
    colors = light[1]['color']
    i = 0
    while i < 3:
        ret[i]=colors[i]*ctheta
        i+=1
    return ret

def gen_ispec(matrix, point, light, normal):
    ret = [0,0,0]
    colors = light[1]['color']
    l = normalize(light[1]['location'])
    n = normalize(normal)
    r = sub(mult(n, 2*dot(n, l)), l)
    v = [0,0,1]
    calpha = dot(r,v)
    i = 0
    while i < 3:
        ret[i]=colors[i]*calpha
        i+=1
    return ret

def sortPoints(matrix, point):
    y1 = matrix[point][1]
    y2 = matrix[point+1][1]
    y3 = matrix[point+2][1]
    if y1==min(y1,y2,y3):
	if y2 < y3:
            return matrix[point],matrix[point+1],matrix[point+2]
        else:
            return matrix[point],matrix[point+2],matrix[point+1]
    elif y2==min(y1, y2, y3):
        if y1 < y3:
            return matrix[point+1],matrix[point],matrix[point+2]
        else:
            return matrix[point+1],matrix[point+2],matrix[point]
    else:
        if y1 < y2:
            return matrix[point+2],matrix[point],matrix[point+1]
        else:
            return matrix[point+2],matrix[point+1],matrix[point]

def scanline_convert(matrix, point, screen, zbuff, color, normal):
    newcolor = gen_color(matrix, point, color, normal)
    a = int(newcolor[0])
    b = int(newcolor[1])
    c = int(newcolor[2])
    coors = sortPoints(matrix, point)
    bx = float(coors[0][0])
    mx = float(coors[1][0])
    tx = float(coors[2][0])
    bz = float(coors[0][2])
    mz = float(coors[1][2])
    tz = float(coors[2][2])
    by = int(coors[0][1])
    my = int(coors[1][1])
    ty = int(coors[2][1])
    y = by
    x0 = bx
    x1 = bx
    z0 = bz
    z1 = bz
    if by==my:
        x1 = mx
        z1 = mz
    while y <= ty:
        x0+= (tx-bx)/(ty-by)
        z0 += (tz-bz)/(ty-by)
        if y < my:
            x1 += (mx-bx)/(my-by)
            z1 += (mz-bz)/(my-by)
        elif y>my:
            x1 += (tx-mx)/(ty-my)
            z1 += (tz-mz)/(ty-my)
        y+=1
        draw_line(int(x0), int(y), int(z0), int(x1), int(y), int(z1), screen, zbuff, [a,b,c])

def add_polygon( polygons, x0, y0, z0, x1, y1, z1, x2, y2, z2 ):
    add_point(polygons, x0, y0, z0);
    add_point(polygons, x1, y1, z1);
    add_point(polygons, x2, y2, z2);

def draw_polygons( matrix, screen, zbuffer, color ):
    if len(matrix) < 2:
        print 'Need at least 3 points to draw'
        return

    point = 0    
    while point < len(matrix) - 2:

        normal = calculate_normal(matrix, point)[:]
        if normal[2] > 0:
            scanline_convert(matrix, point, screen, zbuffer, color, normal)
        point+= 3


def add_box( polygons, x, y, z, width, height, depth ):
    x1 = x + width
    y1 = y - height
    z1 = z - depth

    #front
    add_polygon(polygons, x, y, z, x1, y1, z, x1, y, z);
    add_polygon(polygons, x, y, z, x, y1, z, x1, y1, z);
  
    #back
    add_polygon(polygons, x1, y, z1, x, y1, z1, x, y, z1);
    add_polygon(polygons, x1, y, z1, x1, y1, z1, x, y1, z1);
  
    #right side
    add_polygon(polygons, x1, y, z, x1, y1, z1, x1, y, z1);
    add_polygon(polygons, x1, y, z, x1, y1, z, x1, y1, z1);
    #left side
    add_polygon(polygons, x, y, z1, x, y1, z, x, y, z);
    add_polygon(polygons, x, y, z1, x, y1, z1, x, y1, z);
  
    #top
    add_polygon(polygons, x, y, z1, x1, y, z, x1, y, z1);
    add_polygon(polygons, x, y, z1, x, y, z, x1, y, z);
    #bottom
    add_polygon(polygons, x, y1, z, x1, y1, z1, x1, y1, z);
    add_polygon(polygons, x, y1, z, x, y1, z1, x1, y1, z1);

def add_sphere( edges, cx, cy, cz, r, step ):
    points = generate_sphere(cx, cy, cz, r, step)
    num_steps = int(1/step+0.1)
    
    lat_start = 0
    lat_stop = num_steps
    longt_start = 0
    longt_stop = num_steps

    num_steps+= 1
    for lat in range(lat_start, lat_stop):
        for longt in range(longt_start, longt_stop):
            
            p0 = lat * (num_steps) + longt
            p1 = p0+1
            p2 = (p1+num_steps) % (num_steps * (num_steps-1))
            p3 = (p0+num_steps) % (num_steps * (num_steps-1))

            if longt != num_steps - 2:
	        add_polygon( edges, points[p0][0],
		             points[p0][1],
		             points[p0][2],
		             points[p1][0],
		             points[p1][1],
		             points[p1][2],
		             points[p2][0],
		             points[p2][1],
		             points[p2][2])
            if longt != 0:
	        add_polygon( edges, points[p0][0],
		             points[p0][1],
		             points[p0][2],
		             points[p2][0],
		             points[p2][1],
		             points[p2][2],
		             points[p3][0],
		             points[p3][1],
		             points[p3][2])

def generate_sphere( cx, cy, cz, r, step ):
    points = []
    num_steps = int(1/step+0.1)
    
    rot_start = 0
    rot_stop = num_steps
    circ_start = 0
    circ_stop = num_steps
            
    for rotation in range(rot_start, rot_stop):
        rot = step * rotation
        for circle in range(circ_start, circ_stop+1):
            circ = step * circle

            x = r * math.cos(math.pi * circ) + cx
            y = r * math.sin(math.pi * circ) * math.cos(2*math.pi * rot) + cy
            z = r * math.sin(math.pi * circ) * math.sin(2*math.pi * rot) + cz

            points.append([x, y, z])
            #print 'rotation: %d\tcircle%d'%(rotation, circle)
    return points
        
def add_torus( edges, cx, cy, cz, r0, r1, step ):
    points = generate_torus(cx, cy, cz, r0, r1, step)
    num_steps = int(1/step+0.1)
    
    lat_start = 0
    lat_stop = num_steps
    longt_start = 0
    longt_stop = num_steps
    
    for lat in range(lat_start, lat_stop):
        for longt in range(longt_start, longt_stop):

            p0 = lat * (num_steps) + longt;
            if (longt == num_steps - 1):
	        p1 = p0 - longt;
            else:
	        p1 = p0 + 1;
            p2 = (p1 + num_steps) % (num_steps * num_steps);
            p3 = (p0 + num_steps) % (num_steps * num_steps);

            add_polygon(edges,
                        points[p0][0],
                        points[p0][1],
                        points[p0][2],
                        points[p3][0],
                        points[p3][1],
                        points[p3][2],
                        points[p2][0],
                        points[p2][1],
                        points[p2][2] )
            add_polygon(edges,
                        points[p0][0],
                        points[p0][1],
                        points[p0][2],
                        points[p2][0],
                        points[p2][1],
                        points[p2][2],
                        points[p1][0],
                        points[p1][1],
                        points[p1][2] )

def generate_torus( cx, cy, cz, r0, r1, step ):
    points = []
    num_steps = int(1/step+0.1)
    
    rot_start = 0
    rot_stop = num_steps
    circ_start = 0
    circ_stop = num_steps
    
    for rotation in range(rot_start, rot_stop):
        rot = step * rotation
        for circle in range(circ_start, circ_stop):
            circ = step * circle

            x = math.cos(2*math.pi * rot) * (r0 * math.cos(2*math.pi * circ) + r1) + cx;
            y = r0 * math.sin(2*math.pi * circ) + cy;
            z = -1*math.sin(2*math.pi * rot) * (r0 * math.cos(2*math.pi * circ) + r1) + cz;

            points.append([x, y, z])
    return points

def add_circle( points, cx, cy, cz, r, step ):
    x0 = r + cx
    y0 = cy
    t = step

    while t <= 1.00001:
        x1 = r * math.cos(2*math.pi * t) + cx;
        y1 = r * math.sin(2*math.pi * t) + cy;
        add_edge(points, x0, y0, cz, x1, y1, cz)
        x0 = x1
        y0 = y1
        t+= step

#generates an xy circle to rotate
def xy_circle(r, step):
    points = []
    theta = step
    x=0
    y=0
    z=0
    while theta <= 1.0001:
       x=r*math.cos(2*math.pi * theta)
       y=r*math.sin(2*math.pi * theta)
       theta+=step
       points.append([x,y,z,1])
    return points

def add_curve( points, x0, y0, x1, y1, x2, y2, x3, y3, step, curve_type ):

    xcoefs = generate_curve_coefs(x0, x1, x2, x3, curve_type)[0]
    ycoefs = generate_curve_coefs(y0, y1, y2, y3, curve_type)[0]

    t = step
    while t <= 1.00001:
        x = xcoefs[0] * t*t*t + xcoefs[1] * t*t + xcoefs[2] * t + xcoefs[3]
        y = ycoefs[0] * t*t*t + ycoefs[1] * t*t + ycoefs[2] * t + ycoefs[3]
                
        add_edge(points, x0, y0, 0, x, y, 0)
        x0 = x
        y0 = y
        t+= step

def gen_bezier3( points, x0, y0, z0, x1, y1, z1, x2, y2, z2, x3, y3, z3, r, step, screen, zbuffer, color):
    xcoefs = generate_curve_coefs(x0, x1, x2, x3, 'bezier')[0]
    ycoefs = generate_curve_coefs(y0, y1, y2, y3, 'bezier')[0]
    zcoefs = generate_curve_coefs(z0, z1, z2, z3, 'bezier')[0]
    t = step
    while t <= 1.00001:
        x = xcoefs[0] * t*t*t + xcoefs[1] * t*t + xcoefs[2] * t + xcoefs[3]
        dx = 3*xcoefs[0] * t*t + 2*xcoefs[1] * t + xcoefs[2]
        y = ycoefs[0] * t*t*t + ycoefs[1] * t*t + ycoefs[2] * t + ycoefs[3]
        dy = 3*ycoefs[0] * t*t + 2*ycoefs[1] * t + ycoefs[2]
        z = zcoefs[0] * t*t*t + zcoefs[1] * t*t + zcoefs[2] * t + zcoefs[3]
        dz = 3*zcoefs[0] * t*t + 2*zcoefs[1] * t + zcoefs[2]
        xyz = normalize([x,y,z])
        thetax = math.acos(xyz[0])
        thetay = math.acos(xyz[1])
        thetaz = math.acos(xyz[2])
        circle = xy_circle(r, step)
        matrix_mult(make_rotX(thetax),circle)
        matrix_mult(make_rotY(thetay),circle)
        matrix_mult(make_rotZ(thetaz),circle)
        matrix_mult(make_translate(x,y,z),circle)
        print_matrix(circle)
        i = 1
        while i < len(circle):
            draw_line( circle[i-1][0], circle[i-1][1], circle[i-1][2], circle[i][0], circle[i][1], circle[i][2], screen, zbuffer, color )
            i+=1
        draw_line( circle[-1][0], circle[-1][1], circle[-1][2], circle[0][0], circle[0][1], circle[0][2], screen, zbuffer, color )
        t+=step
    return points

def draw_lines( matrix, screen, zbuffer, color ):
    if len(matrix) < 2:
        print 'Need at least 2 points to draw'
        return
    
    point = 0
    while point < len(matrix) - 1:
        draw_line( int(matrix[point][0]),
                   int(matrix[point][1]),
                   matrix[point][2],
                   int(matrix[point+1][0]),
                   int(matrix[point+1][1]),
                   matrix[point+1][2],
                   screen, zbuffer, color)    
        point+= 2
        
def add_edge( matrix, x0, y0, z0, x1, y1, z1 ):
    add_point(matrix, x0, y0, z0)
    add_point(matrix, x1, y1, z1)
    
def add_point( matrix, x, y, z=0 ):
    matrix.append( [x, y, z, 1] )
    



def draw_line( x0, y0, z0, x1, y1, z1, screen, zbuffer, color ):
    #swap points if going right -> left
    if x0 > x1:
        xt = x0
        yt = y0
        zt = z0
        x0 = x1
        y0 = y1
        z0 = z1
        x1 = xt
        y1 = yt
        z1 = zt

    x = x0
    y = y0
    z = z0
    A = 2 * (y1 - y0)
    B = -2 * (x1 - x0)
    wide = False
    tall = False

    if ( abs(x1-x0) >= abs(y1 - y0) ): #octants 1/8
        wide = True
        loop_start = x
        loop_end = x1
        dx_east = dx_northeast = 1
        dy_east = 0
        d_east = A
        distance = x1 - x
        if ( A > 0 ): #octant 1
            d = A + B/2
            dy_northeast = 1
            d_northeast = A + B
        else: #octant 8
            d = A - B/2
            dy_northeast = -1
            d_northeast = A - B

    else: #octants 2/7
        tall = True
        dx_east = 0
        dx_northeast = 1
        distance = abs(y1 - y)
        if ( A > 0 ): #octant 2
            d = A/2 + B
            dy_east = dy_northeast = 1
            d_northeast = A + B
            d_east = B
            loop_start = y
            loop_end = y1
        else: #octant 7
            d = A/2 - B
            dy_east = dy_northeast = -1
            d_northeast = A - B
            d_east = -1 * B
            loop_start = y1
            loop_end = y

    while ( loop_start < loop_end ):
        plot( screen, zbuffer, color, x, y, z )
        if ( (wide and ((A > 0 and d > 0) or (A < 0 and d < 0))) or
             (tall and ((A > 0 and d < 0) or (A < 0 and d > 0 )))):
            x+= dx_northeast
            y+= dy_northeast
            d+= d_northeast
        else:
            x+= dx_east
            y+= dy_east
            d+= d_east
        loop_start+= 1

    plot( screen, zbuffer, color, x, y, z )

    
