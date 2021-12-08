import math

def RGB2XYZ(tup0):
    (R, G, B) = tup0
    # R, G and B (Standard RGB) input range=0 to 255

    red = R / 255  # normalize red
    green = G / 255  # normalize green
    blue = B / 255  # normalize blue

    if red > 0.04045:
        red_linear = math.pow(((red + 0.055) / 1.055), 2.4)
    else:
        red_linear = red / 12.92
    if green > 0.04045:
        green_linear = math.pow(((green + 0.055) / 1.055), 2.4)
    else:
        green_linear = green / 12.92
    if blue > 0.04045:
        blue_linear = math.pow(((blue + 0.055) / 1.055), 2.4)
    else:
        blue_linear = blue / 12.92
    # these red, blue, and green are linearized r,g,b. They are gamma adjusted, where gamma is 2.4 per sRGB

    red_linear = red_linear * 100
    green_linear = green_linear * 100
    blue_linear = blue_linear * 100

    X = red_linear * 0.4124 + green_linear * 0.3576 + blue_linear * 0.1805
    Y = red_linear * 0.2126 + green_linear * 0.7152 + blue_linear * 0.0722
    Z = red_linear * 0.0193 + green_linear * 0.1192 + blue_linear * 0.9505

    # Reference-X, Y and Z refer to specific illuminants and observers.
    # Use D65 Reference X, Y,Z values for 2 degrees (CIE 1931) for sRGB
    # Reference_X=95.047 #D65 Reference X, Y,Z values for 2 degrees (CIE 1931)
    # Reference_Y=100.000 #D65 Reference X, Y,Z values for 2 degrees (CIE 1931)
    # Reference_Z=108.883 #D65 Reference X, Y,Z values for 2 degrees (CIE 1931)
    Reference_X = 100
    Reference_Y = 100
    Reference_Z = 100

    X = X / Reference_X
    Y = Y / Reference_Y
    Z = Z / Reference_Z

    return (X, Y, Z)
def XYZ2RGB(tup0):
    (X, Y, Z) = tup0

    red_linear = X * 3.2404542 + (-1.5371385) * Y + (-0.4985314) * Z
    green_linear = X * (-0.9692660) + Y * 1.8760108 + Z * 0.0415560
    blue_linear = X * 0.0556434 + Y * (-0.2040259) + Z * 1.0572252

    if red_linear <= 0: red_linear = 0
    if green_linear <= 0: green_linear = 0
    if blue_linear <= 0: blue_linear = 0

    if (red_linear > 0.0031308):
        red = 1.055 * math.pow(red_linear, 0.4166666) - 0.055
    elif red_linear >= 0:
        red = 12.92 * red_linear
    else:
        red = 0

    if (green_linear > 0.0031308):
        green = 1.055 * math.pow(green_linear, 0.4166666) - 0.055
    elif green_linear >= 0:
        green = 12.92 * green_linear
    else:
        green = 0

    if (blue_linear > 0.0031308):
        blue = 1.055 * math.pow(blue_linear, 0.4166666) - 0.055
    elif blue_linear >= 0:
        blue = 12.92 * blue_linear
    else:
        blue = 0

    if red >= 1: red = 1
    if green >= 1: green = 1
    if blue >= 1: blue = 1

    R = int(red * 255)
    G = int(green * 255)
    B = int(blue * 255)

    return (R, G, B)

def colormix_xyz(color1,color2,r1,r2):
#To mix colors one has to use linear color space like XYZ or linearized rgb (adjusted for gamma of the sRGB color space).
#This function uses the linear XYZ interpolation

    perc1=r1/(r1+r2)
    perc2=1-perc1
    outputimage=Image.new('RGB',(2*image_width,image_height))
    outputimage_pixelmap=outputimage.load()
    col1=RGB2XYZ(color1)
    col2=RGB2XYZ(color2)
    xyz_linear=((col1[0]*perc1+col2[0]*perc2),(col1[1]*perc1+col2[1]*perc2),(col1[2]*perc1+col2[2]*perc2))
    rgb_linearmix=XYZ2RGB(xyz_linear)
    #Next 2 nested for loops output 2 images side by side. First image is linear optical mix of color1 and color2. Second image is the resultant rgb_linearmix
    for j in range(image_height): #randomize the distribution of color1 and color2 proportional to r1 and r2
        for i in range(image_width):
            rn=random.random()
            if rn<=perc1:
                pixcolor=color1
            else:
                pixcolor=color2
            outputimage_pixelmap[i,j]=pixcolor
    for j in range(image_height):
        for i in range(image_width+1,2*image_width):
            outputimage_pixelmap[i,j]=rgb_linearmix
    outputimage.save('temp-xyzavg.png')
    outputimage.show()
