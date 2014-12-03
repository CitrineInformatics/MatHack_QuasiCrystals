import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import host_subplot
import mpl_toolkits.axisartist as AA
import numpy as np
import math
from numpy.linalg import inv
from pyhull.delaunay import DelaunayTri
from matplotlib.path import Path
import matplotlib.patches as patches
from pymatgen.command_line.gulp_caller import *


pi=math.pi

'''
Materials Hackathon
Fall MRS 2014

Wenhao Sun
MIT
'''



def dc2rs(DC,N,bonddist=1):
    x=0.0
    y=0.0
    for nn in range(0,N):
        x=x+bonddist*DC[nn]*np.cos(nn*2*pi/N)
        y=y+bonddist*DC[nn]*np.sin(nn*2*pi/N)
    return [x,y]

def getdualcoord(N,nstrips,C,joggle):
    DualCoord=np.zeros(N)
    for nn in range(0,N):
        if N%2==0:
            angle=2*pi*nn/N+0.03
        else:
            angle=2*pi*nn/N
        m=np.tan(angle)
        ## Generate all intersects
        if angle==0:
            db=1.00
        elif 0<angle<=pi/2 or 3*pi/2<angle<2*pi:
            db=1/np.sin(pi/2-angle)
        else:
            db=1/np.sin(angle-pi/2)
        minb=-nstrips*db/2+joggle[nn]
        #Find which strip it lies in
        dy=C[1]-m*C[0]-minb
        strip=np.ceil(dy/db)
        if dy/db<0:
            np.floor(dy/db)
        if pi/2<angle<=3*pi/2:
            strip=nstrips-(strip-1)
        else:
            pass
        #print m,minb,strip
        DualCoord[nn]=strip
    #print DualCoord
    return DualCoord

#Initialize, make size and dimensionality variables
def quasi(N, nstrips=10, disp=False,tojoggle=False):
    #Can add a joggle to get a more stochastic quasicrystal 
    if tojoggle==False:
        joggle=np.zeros(N)
    else:
        joggle=np.random.rand(N)
    
    #N = characteristic dimensionality
    if nstrips%2==0:
        nstrips+=1
    #Find all lines
    setoflines=[] #Parallel lines
    for nn in range(0,N):
        ## Generate slopes
        lines=[]
        if N%2==0:
            angle=2*pi*nn/N+0.03
        else:
            angle=2*pi*nn/N
        m=np.tan(angle)
        ## Generate all intersects
        if 0<angle<pi/2 or 3*pi/2<angle<2*pi:
            db=1/np.sin(pi/2-angle)
        else:
            db=1/np.sin(angle-pi/2)
        minb=-nstrips*db/2+joggle[nn]
        for bb in range(0,nstrips+1):
            b=minb+db*bb
            line=[m,b]
            lines.append(line)
        setoflines.append(lines)

    """
    All lines
    xx = np.arange(-5, 5, 0.01)
    for lines in setoflines:
        for line in lines:
            [m,b]=line
            plt.plot(xx,m*xx+b)
            plt.ylim(-5,5)
    plt.draw()
    plt.show()
    """
    #Find all intersects of lines
    intersects=[]
    for n1 in range(0,N-1):
        for n2 in range(n1+1,N):
            for line1 in setoflines[n1]:
                for line2 in setoflines[n2]:
                    A=np.array([[-line1[0],1],[-line2[0],1]])
                    B=np.array([[line1[1]],[line2[1]]])
                    I=np.dot(inv(A),B)
                    intersects.append([float(I[0]),float(I[1])])

    NewIntersects=[]
    #Keep only intersects inside a radius
    r0=np.floor(nstrips/2)
    for I in intersects:
        r=np.sqrt(I[0]**2+I[1]**2)
        if r<r0:
            NewIntersects.append(I)
    intersects=NewIntersects
    
    
    if disp==True:
        xx = np.arange(-10, 10, 0.01)
        scatterx=[]
        scattery=[]
        for I in intersects:
            scatterx.append(I[0])
            scattery.append(I[1])
        plt.scatter(scatterx,scattery)
        for lines in setoflines:
            for line in lines:
                [m,b]=line
                plt.plot(xx,m*xx+b)
                plt.axis([-r0-2,r0+2,-r0-2,r0+2])
                plt.draw()
        plt.show()
        
    #Make Delaunay Triangulation. 
    #Ideally you would get one point per enclosed polgon, but I actually don't know how to code this
    #This way you get multiple points per polygon, which I will then sort. 
    #It's inefficient but who cares it's a hackathon. 
    tri = DelaunayTri(intersects)
    
    #Get convex centroids
    centroids=[]
    for points in tri.vertices:
        center=np.zeros(2)
        for p in points:
            center[0]+=intersects[p][0]
            center[1]+=intersects[p][1]
        centroid=(center[0]/3,center[1]/3)
        centroids.append(centroid)
    #print len(centroids)
    
    #Keep only centroids inside a radius
    NewCentroids=[]
    r0=r0-0.5
    for I in centroids:
        r=np.sqrt(I[0]**2+I[1]**2)
        if r<r0:
            NewCentroids.append(I)
    centroids=NewCentroids
    if disp==True:
        xx = np.arange(-r0-2, r0+2, 0.01)
        scatterx=[]
        scattery=[]
        for I in centroids:
            scatterx.append(I[0])
            scattery.append(I[1])
        plt.scatter(scatterx,scattery)
        for lines in setoflines:
            for line in lines:
                [m,b]=line
                plt.plot(xx,m*xx+b)
                plt.axis([-r0-2,r0+2,-r0-2,r0+2])
                plt.draw()
        plt.show()
        
    #Get N-dimensional dual-space 'coordinates' of centroids
    DualCoords=[]
    UniqueCentroids=[]
    r0=r0-0.5
    for C in centroids:
        DualCoord=getdualcoord(N,nstrips,C,joggle)
        #if any(DualCoord[x]>nstrips for x in range(0,N)) or any(DualCoord[x]<0 for x in range(0,N)):
            #print "Discard"
        #    continue
        """
        FOR DEBUGGING
        
        lines=[]
        for bb in range(0,nstrips+1):
            b=minb+db*bb
            line=[m,b]
            lines.append(line)
        plt.scatter(C[0],C[1])
        xx = np.arange(-10, 10, 0.01)
        for lines in setoflines:
            for line in lines:
                [m,b]=line
                plt.plot(xx,m*xx+b)
                plt.ylim(-10,10)
                plt.draw()
        plt.show()
        """
        #Make sure it's not doublecounted
        r=np.sqrt(C[0]**2+C[1]**2)
            
        if not DualCoord==None and not any((DualCoord == x).all() for x in DualCoords) and r<r0:
           DualCoords.append(DualCoord)
           UniqueCentroids.append(C) 
    
    if disp==True:
        xx = np.arange(-5, 5, 0.01)
        scatterx=[]
        scattery=[]
        for I in UniqueCentroids:
            scatterx.append(I[0])
            scattery.append(I[1])
        plt.scatter(scatterx,scattery)
        for lines in setoflines:
            for line in lines:
                [m,b]=line
                plt.plot(xx,m*xx+b)
                plt.ylim(-5,5)
                plt.draw()
        plt.show()
    #print len(DualCoords)    

    #Do dual space to real space transformation
    RealCoords=[]
    for DC in DualCoords:
        rs=dc2rs(DC,N)
        #print DC, rs
        RealCoords.append(rs)
    #print len(RealCoords)
    
    #Get connecting lines - draw a few extra on the outside
    intersects=[]
    tilelines=[]
    fig = plt.figure()
    for n1 in range(0,N-1):
        for n2 in range(n1+1,N):
            for line1 in setoflines[n1]:
                for line2 in setoflines[n2]:
                    A=np.array([[-line1[0],1],[-line2[0],1]])
                    B=np.array([[line1[1]],[line2[1]]])
                    I=np.dot(inv(A),B)
                    newI=[float(I[0]),float(I[1])]
                    if newI in NewIntersects:
                        avgslope=(line1[0]+line2[0])/2
                        perpslope=-1/avgslope
                        fourpoints=[]
                        for jj in [avgslope, perpslope]:
                            for ii in [-1,1]:
                                t=np.arctan(jj)
                                x=I[0]+0.05*np.cos(t)*ii
                                y=I[1]+0.05*np.sin(t)*ii
                                fourpoints.append([float(x),float(y)])
    
                        rfp=[]
                        for p in fourpoints:
                            dc=getdualcoord(N,nstrips,p,joggle) 
                            rs=dc2rs(dc,N) 
                            rfp.append(rs)
                        tile=[]
                        tile.append(rfp[0])
                        tile.append(rfp[2])
                        tile.append(rfp[1])
                        tile.append(rfp[3])
                        tile.append(rfp[0])
                        codes = [Path.MOVETO,
                         Path.LINETO,
                         Path.LINETO,
                         Path.LINETO,
                         Path.CLOSEPOLY,
                        ]
                        path = Path(tile, codes)
                        ax = fig.add_subplot(111)
                        patch = patches.PathPatch(path, facecolor='none', lw=0.5)
                        ax.add_patch(patch)
    #Construct Real Space Quasicrystal
    #disp=True
    if disp==True:
        scatterx=[]
        scattery=[]
        for I in RealCoords:
            scatterx.append(I[0])
            scattery.append(I[1])
        plt.scatter(scatterx,scattery,s=40)
        plt.axis([-nstrips-2,nstrips+2,-nstrips-2,nstrips+2])
        plt.show()
    ### I think it is drawing some extra lines because the 'strip' (the Dual coordinate) of very edge atoms are not well accounted for. 
    return RealCoords

def get_quasicrystal_energy(coordinates,scale=5):
    #Make box very large (100Ax100A), so that there are no interactions across the periodic boundary
    #Using GULP to calculate energies
    #Caller from pymatgen.command_line.gulp_caller
    #All atoms set to be silver, because there's an example potential on the GULP website for it. Would best have been aluminum probably
    #Create surface unit cell input
    
    gin="single\n"
    gin+="scell \n"
    S=str(scale*100)
    gin+=S+" "+S+" 90 \n"
    gin+="sfractional region 1 \n"
    for c in coordinates:
        gin+="Ag  core   "+str(c[0]/100)+"  "+str(c[1]/100)+"  0 0 1 0 \n" 
    gin+="""species
Ag core  0.000
morse 
Ag core Ag core 0.6721 1.8260 2.5700 0.0 5.542
cutp 5.542 voter 
manybody
Ag core Ag core  0.0 5.5420
eam_density voter
Ag core   1.0 3.9060
eam_potential_shift 
Ag core Ag core 10.4262977 3.9060 5.542
eam_functional numeric
Ag core ag.dbout"""
    gc = GulpCaller()
    gout = gc.run(gin)
    for line in gout.split("\n"):
        if "Total lattice energy" in line and "eV" in line:
            energy=line.split()
    return float(energy[4])

def energy_conjugate_gradient(coordinates): 
    #Change scale until lowest energy converged
    # 6 AM ... pretty tired
    minE=100
    minS=0
    for scale in range(1,20):
        s=3+0.2*scale
        E=get_quasicrystal_energy(coordinates,s)
        if E<minE:
            minE=E
            minS=s
    #print minS,minE
    x0=minS #scale
    TOL=0.001*len(coordinates)
    #print TOL
    n=1
    
    def f_(coordinates,x0):
        A=get_quasicrystal_energy(coordinates,x0+0.001)
        B=get_quasicrystal_energy(coordinates,x0)
        return (A-B)/0.001
        
    while n<=100:
        x1 = x0 - f_(coordinates,x0)*0.001
        E=get_quasicrystal_energy(coordinates,x0)
        #print x0,E
        if np.abs(x1 - x0) < TOL:
            return x1, E
        else:
            x0 = x1
    return False


     
if __name__ == '__main__':
    #quasi(N,nstrips)
    for dim in [3,5,7,8,9,10]:
        print dim
        for nstrips in [13,15]:
            rc=quasi(dim,nstrips,disp=False,tojoggle=False)
            [s,E]=energy_conjugate_gradient(rc)
            print nstrips,len(rc),E
    
"""
    #Algorithm Figures
    quasi(5,9,disp=True)
    
    #Dimensionality Figures
    for dim in [3,5,7,8,9,10]:
        print dim
        for nstrips in [9]:
            rc=quasi(dim,nstrips,disp=False)#Need to set the last figure for disp=true
    
    #Demonstrate joggling
    quasi(5,9,disp=True,tojoggle=True)
    
    ##Get all data
    for dim in [3,5,7,8,9,10]:
        print dim
        for nstrips in [5,7,11]:
            rc=quasi(dim,nstrips,disp=False,tojoggle=False)
            [s,E]=energy_conjugate_gradient(rc)
            print nstrips,len(rc),E
    """