from config import *

def create_new_shape_gene():
    genestring='' #genestring is a bunch of 1s and 0s that correspond to (cx,cy,col_index) which makes one shape
    for i in range(no_of_shapes):
        b_cx = int(random.random() * 2**pos_bit_depth)
        b_cy = int(random.random() * 2**pos_bit_depth)
        zz=int(random.random()*2**pos_bit_depth)
        col_index=int(random.random()*n_color)
        if col_index==n_color:col_index=int(n_color-1)
        gene_cx = bin(b_cx)[2:].zfill(pos_bit_depth) #convert decimal to binary string, eliminate the first 2 letters, and fill the beginning with zeros
        gene_cy = bin(b_cy)[2:].zfill(pos_bit_depth) #convert decimal to binary string, eliminate the first 2 letters, and fill the beginning with zeros
        gene_zz=bin(zz)[2:].zfill(pos_bit_depth)
        gene_col=bin(col_index)[2:].zfill(col_bit_depth) #convert decimal to binary string, eliminate the first 2 letters, and fill the beginning with zeros
        indiv_shape=gene_cx+gene_cy+gene_zz+gene_col #this is a string that is xx characters long, looks like a binary number
        genestring=genestring+indiv_shape
    return genestring #genestring is a bunch of 1s and 0s

def create_new_shape_gene_v2():
    gene_shape_list = []  # genestring is a bunch of 1s and 0s that correspond to (cx,cy,col_index) which makes one shape
    for i in range(no_of_shapes):
        nsg=create_new_shape_gene()
        dgf=decode_gene_to_float(nsg)
        gene_shape_list=gene_shape_list+list(dgf)
    return gene_shape_list #this is in format [cx0,cy1,zz0,col_index0,cx1,cy1,zz1,col_index1,....) cx and cy values
    # are floats between 0 and 1

def create_new_binary_population(n_pop):
    pop_list=[]
    for i in range(n_pop):
        pop_list.append(create_new_shape_gene())
    return pop_list #poplist is in format [ '1010001010...0101', '10110011...01011',....]

def create_new_population_v2(n_pop):
    pop_list = []
    for i in range(n_pop):
        pop_list.append(create_new_shape_gene_v2())
    return pop_list

def DNAmutate(bin_genestring,mutation_prob):

    def mutate(gen,p0):
        prob=random.random()
        if prob<=p0:
            return str((int(gen)+1)%2)
        else:
            return gen

    n=len(bin_genestring)
    mutation_ct=int(n/20) #200 is an arbitrary number
    mutation_index_list=[]
    for i in range(mutation_ct):
        mutation_index_list.append(int(random.random()*n))
    ci_list=sorted(mutation_index_list)
    genestring=bin_genestring
    for index in ci_list:
        newgene=mutate(bin_genestring[index],mutation_prob)
        new_genestring=genestring[:index]+newgene+genestring[index+1:]
        genestring=new_genestring
    return genestring

def DNAexchange(indiv1,indiv2):
    max_n=min(len(indiv1),len(indiv2))
    mem1=indiv1
    mem2=indiv2
    crossover_ct=int(max_n/800) #20 is an arbitrary number
    if crossover_ct==0: crossover_ct=1
    crossover_index_list=[]
    for i in range(crossover_ct):
        crossover_index_list.append(int(random.random()*max_n))
    ci_list=sorted(crossover_index_list)
    for index in ci_list:
        new_mem1=mem1[:index]+mem2[index:]
        new_mem2=mem2[:index]+mem1[index:]
        mem1=new_mem1
        mem2=new_mem2
    return (mem1,mem2)

def evolve_binary_population(bin_population):
    p_mutate=0.05
    n_pop=len(bin_population)
    n_mates=int(n_pop/2)
    newpopulation=[]
    for i in range(n_mates):
        indiv_1=bin_population[2*i]
        indiv_2=bin_population[2*i+1]
        offsprings=DNAexchange(indiv_1,indiv_2)
        newpopulation.append(DNAmutate(offsprings[0],p_mutate)) #add a mutated version of the new offspring
        newpopulation.append(DNAmutate(offsprings[1], p_mutate))
    if n_pop%2==1:
        additional_offsprings=DNAexchange(bin_population[-1],newpopulation[-1])
        newpopulation.pop(-1) #this removes the last member of the new population
        newpopulation.append(DNAmutate(additional_offsprings[0], p_mutate))  # add a mutated version of the new offspring
        newpopulation.append(DNAmutate(additional_offsprings[1], p_mutate))
    return newpopulation

def evolve_nb_population(nb_population):
    n_pop=len(nb_population)
    n_mates=int(n_pop/2)
    newpopulation=[]
    for i in range(n_mates):
        indiv_1=nb_population[2*i]
        indiv_2=nb_population[2*i+1]
        offsprings=DNAexchange(indiv_1,indiv_2)
        newpopulation.append(offsprings[0])
        newpopulation.append(offsprings[1])
    if n_pop%2==1:
        additional_offsprings=DNAexchange(bin_population[-1],newpopulation[-1])
        newpopulation.pop(-1) #this removes the last member of the new population
        newpopulation.append(additional_offsprings[0])  # add a mutated version of the new offspring
        newpopulation.append(additional_offsprings[1])
    return newpopulation

def decode_gene_to_float(gene):
    cx = int('0b' + gene[0:pos_bit_depth], 2) / (2**pos_bit_depth)  # this converts cx to range [0,1]
    cy = int('0b' + gene[pos_bit_depth:pos_bit_depth * 2], 2) / (2**pos_bit_depth)  # this converts cy to range [0,1]
    zz=int('0b'+gene[pos_bit_depth * 2:pos_bit_depth * 3],2)
    col_index = int('0b' + gene[-col_bit_depth:], 2) % n_color # in case of mutations in genes, the index can not exceed total number of colors
    return (cx,cy,zz,col_index)

def decode_gene_seq(glist,w,h,x0,y0):
    #glist is a string consisting of 1s and 0s. It is an encryption of shapes of n_shape quantity.
    #w, h, x0, y0 are all in inches or mm with respect to the zone upper left corner
    n_string=3*pos_bit_depth+col_bit_depth #this is the length of a genestring for each individual
    n_gene=len(glist)
    n_shape=n_gene/n_string
    shape_list=[]
    if int(n_shape)!= n_shape:
        print ('****DECODING ERROR****')
        return
    for nn in range(int(n_shape)):
        i=int(nn*n_string)
        indiv_shape=glist[i:i+n_string]
        decoded_f=decode_gene_to_float(indiv_shape)
        (cx,cy,zz,col_index)=(decoded_f[0]*w+x0,decoded_f[1]*h+y0,decoded_f[2],decoded_f[3]) #decoded shape is within the zone
        #but cx and cy are in inches or mm in global coordinates.
        #print ('cx,cy: ', cx, ' , ',cy)
        new_brush=brushstroke((cx,cy),zz,pigments[col_index])
        shape_list.append(new_brush)
    return shape_list #shape_list is a collection of brush objects.

def decode_gene_seq_v2(gseq,w,h,x0,y0):
    #gseq is a list that looks like  [cx0,cy1,zz0,col_index0,cx1,cy1,zz1,col_index1,....).
    #It codes for shapes of n_shape quantity.
    #w, h, x0, y0 are all in inches or mm with respect to the zone upper left corner
    n_string=4 #this is the length of a genestring for each individual
    n_gene=len(gseq)
    n_shape=n_gene/n_string
    shape_list=[]
    if int(n_shape)!= n_shape:
        print ('****DECODING ERROR****')
        return
    for nn in range(int(n_shape)):
        i=int(nn*n_string)
        decoded_f=gseq[i:i+n_string]
        (cx,cy,zz,col_index)=(decoded_f[0]*w+x0,decoded_f[1]*h+y0,decoded_f[2],decoded_f[3]) #decoded shape is within the zone
        #but cx and cy are in inches or mm in global coordinates.
        #print ('cx,cy: ', cx, ' , ',cy)
        new_brush=brushstroke((cx,cy),zz,pigments[col_index])
        shape_list.append(new_brush)
    return shape_list #shape_list is a collection of brush objects.

def eval_zone(binarystring,zone):
    ((x1, y1), (x2, y2)) = zone.coord
    zone_width = x2 - x1  # inches or mm
    zone_height = y2 - y1  # inches or mm
    ((xi1, yi1), (xi2, yi2)) = zone.pix_coord
    zone_width_pix=xi2-xi1
    zone_height_pix=yi2-yi1
    new_shape_list = decode_gene_seq_v2(binarystring, zone_width, zone_height, x1, y1)
    updated_zone_shape_list = zone.canvas_shape_list + new_shape_list
    target=np.copy(zone.target_image)
    newcanvas=np.copy(zone.canvas_image)
    slist=[]
    for shapes in updated_zone_shape_list:
        z=shapes.z
        pixlist=shapes.exportpixels()
        slist.append((z,pixlist))
    pixprioritylist=sorted(slist,key=lambda x:x[0])
    for z_sorted_pixels in pixprioritylist:
        brushlayer=z_sorted_pixels[1]
        for newpixel in brushlayer:
            i=newpixel[0]-xi1
            j=newpixel[1]-yi1
            if (i>=0 and i<zone_width_pix) and (j>=0 and j<zone_height_pix):
                bg_col =newcanvas[j][i]
                fg_col=newpixel[2]
                new_color=alpha_over(fg_col,bg_col)
                newcanvas[j][i]=new_color[0:3]
    error=roi_errorcalc(target,newcanvas,2)
    zone_error=error[0]+0.1*error[1]
    return zone_error

def eval_zone_v2(geneseq,zone):
    ((x1, y1), (x2, y2)) = zone.coord
    zone_width = x2 - x1  # inches or mm
    zone_height = y2 - y1  # inches or mm
    ((xi1, yi1), (xi2, yi2)) = zone.pix_coord
    zone_width_pix=xi2-xi1
    zone_height_pix=yi2-yi1
    new_shape_list = decode_gene_seq_v2(geneseq, zone_width, zone_height, x1, y1)
    updated_zone_shape_list = zone.canvas_shape_list + new_shape_list
    target=np.copy(zone.target_image)
    newcanvas=np.copy(zone.canvas_image)
    slist=[]
    for shapes in updated_zone_shape_list:
        z=shapes.z
        pixlist=shapes.exportpixels()
        slist.append((z,pixlist))
    pixprioritylist=sorted(slist,key=lambda x:x[0])
    for z_sorted_pixels in pixprioritylist:
        brushlayer=z_sorted_pixels[1]
        for newpixel in brushlayer:
            i=newpixel[0]-xi1
            j=newpixel[1]-yi1
            if (i>=0 and i<zone_width_pix) and (j>=0 and j<zone_height_pix):
                bg_col =newcanvas[j][i]
                fg_col=newpixel[2]
                new_color=alpha_over(fg_col,bg_col)
                newcanvas[j][i]=new_color[0:3]
    error=roi_errorcalc(target,newcanvas,2)
    zone_error=error[0]+0.1*error[1]
    return zone_error

def export_binarystring_image(binarystring,zone):
    ((x1, y1), (x2, y2)) = zone.coord
    zone_width = x2 - x1  # inches or mm
    zone_height = y2 - y1  # inches or mm
    ((xi1, yi1), (xi2, yi2)) = zone.pix_coord
    zone_width_pix = xi2 - xi1
    zone_height_pix = yi2 - yi1
    new_shape_list = decode_gene_seq(binarystring, zone_width, zone_height, x1, y1)
    updated_zone_shape_list = zone.canvas_shape_list + new_shape_list
    target = np.copy(zone.target_image)
    newcanvas = np.copy(zone.canvas_image)
    slist = []
    for shapes in updated_zone_shape_list:
        z = shapes.z
        pixlist = shapes.exportpixels()
        slist.append((z, pixlist))
    pixprioritylist = sorted(slist, key=lambda x: x[0])
    for z_sorted_pixels in pixprioritylist:
        brushlayer = z_sorted_pixels[1]
        for newpixel in brushlayer:
            i = newpixel[0] - xi1
            j = newpixel[1] - yi1
            if (i >= 0 and i < zone_width_pix) and (j >= 0 and j < zone_height_pix):
                bg_col = newcanvas[j][i]
                fg_col = newpixel[2]
                new_color = alpha_over(fg_col, bg_col)
                newcanvas[j][i] = new_color[0:3]
    tempimage = cv2.cvtColor(newcanvas, cv2.COLOR_RGB2BGR)
    filename = str(int(time.time() * 1000)) + '.png'
    cv2.imwrite(filename, tempimage)

def export_geneseq_image(geneseq,zone,filename):
    ((x1, y1), (x2, y2)) = zone.coord
    zone_width = x2 - x1  # inches or mm
    zone_height = y2 - y1  # inches or mm
    ((xi1, yi1), (xi2, yi2)) = zone.pix_coord
    zone_width_pix = xi2 - xi1
    zone_height_pix = yi2 - yi1
    new_shape_list = decode_gene_seq_v2(geneseq, zone_width, zone_height, x1, y1)
    updated_zone_shape_list = zone.canvas_shape_list + new_shape_list
    newcanvas = np.copy(zone.canvas_image)
    slist = []
    for shapes in updated_zone_shape_list:
        z = shapes.z
        pixlist = shapes.exportpixels()
        slist.append((z, pixlist))
    pixprioritylist = sorted(slist, key=lambda x: x[0])
    for z_sorted_pixels in pixprioritylist:
        brushlayer = z_sorted_pixels[1]
        for newpixel in brushlayer:
            i = newpixel[0] - xi1
            j = newpixel[1] - yi1
            if (i >= 0 and i < zone_width_pix) and (j >= 0 and j < zone_height_pix):
                bg_col = newcanvas[j][i]
                fg_col = newpixel[2]
                new_color = alpha_over(fg_col, bg_col)
                newcanvas[j][i] = new_color[0:3]
    tempimage = cv2.cvtColor(newcanvas, cv2.COLOR_RGB2BGR)
    cv2.imwrite(filename, tempimage)
