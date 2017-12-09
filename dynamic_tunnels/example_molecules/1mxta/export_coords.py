# Copyright (c) 2004 Robert L. Campbell
from pymol import cmd
from pymol import stored
import pymol
import glob

def export_coords(object,filename):
    m = cmd.get_model(object)
    util.cbag(object)
    print "Num of atoms: ",len(m.atom)
    coords = []
    vdws = []
    indices = []
    colors = []
    names = []
    resn = []
    resi = []
    symbols = []
    chains = []

    maxidx = -1
    for a in m.atom:
        if a.index > maxidx:
            maxidx = a.index
#   print "Maxidx is ",maxidx

    for a in range(1,maxidx+2):
        coords.append([0,0,0])
        vdws.append(-1) 
        indices.append(-1)
        colors.append([0,0,0])
        resn.append("")
        resi.append(-1)
        names.append("")
        symbols.append("")
        chains.append("")
    
#    print "Created coords:",len(coords)
#    print "Created vdws",len(vdws)
#    print "Created indices",len(indices)
#    print "Created colors",len(colors)
    i = 0
    for a in m.atom:
#        print "Atom:", a
        coords[a.index] = a.coord
        vdws[a.index] = a.vdw
        indices[i] = a.index
        names[i] = a.name
        resi[i] = a.resi
        resn[i] = a.resn
        symbols[i] = a.symbol
        chains[i] = a.chain
        i+=1
#        print "Atom coordinate: ",a.coord, "chain: ", a.chain
#        print "A[",a.index,"],vdw=",a.vdw,"name=",a.name,"resn=",a.resn,"resi=",a.resi
#        print dir(a)

    pymol.stored.scolors = [];
    cmd.iterate(object,"stored.scolors.append([index,color])");
    for c in pymol.stored.scolors:
        t = cmd.get_color_tuple(c[1])
#        print "Color of index ",c[0]," is ",c[1]," that is ",t
        colors[c[0]] = t

    f = open(filename,"w")
    f.write("#columns are:\n")
    f.write("#x y z radiusVDW atomSymbol elemName resName resIdx atomIDX chain\n")
    i = 0
    for a in indices:
        if (a >=0):
#            print i,a,coords[a],vdws[a],colors[a],names[i], "sym:", type(symbols[i])
            f.write("%f %f %f %f %s %s %s %s %d %s\n" % (coords[a][0],coords[a][1],coords[a][2],vdws[a],symbols[i],names[i], resn[i],resi[i],indices[i], chains[i]));
        i+=1
    f.close()

def export_all_pdbs(mask,prefix):
    """
    export_all_pdbs <filemask>,<prefix>

    """
    filelist = glob.glob(mask)
    if filelist:
        filelist.sort()
        for name in filelist:
            cmd.load(name,"tmp")
            export_coords("tmp","%s.%s.coords" % (prefix, name))
            cmd.delete("all")
            print "Saved ",name

cmd.extend('export_coords',export_coords)
cmd.extend('export_all_pdbs',export_all_pdbs)

