import numpy as np
import matplotlib.pyplot as pl

def all_block(tomo,Cmm=None,subdiv=1):
    """return a decomposition of the Cmm and an array of the blocks size

    each block is divided in Cxx,Cyy,Cxy.
    if required those blocks are themselves (regularly) divided
    """

    Nw=tomo.get_Nw()
    Nsubap=tomo.get_Nsubap()


    if(Cmm is None):
        nMeas  =tomo.nMeas()
        Cmm=np.zeros((nMeas,nMeas))
        tomo.matcov(Cmm,0,0,Cmm.shape[1],1)

    blocks=np.zeros((Nw*subdiv*2,Nw*subdiv*2),dtype=np.object)
    size=np.zeros(Nw*subdiv*2,dtype=np.int)


    #for all blocks
    for i in range(Nw):
        x1=np.sum(Nsubap[:i])*2
        x2=x1+Nsubap[i]
        sizeX=Nsubap[i]
        dx=sizeX/subdiv
        if(subdiv*dx<sizeX):
            print "not an integer division: ",sizeX,"/",subdiv
        idx=i*2*subdiv
        idx2=(i*2+1)*subdiv

        size[idx:idx+2*subdiv]=dx

        for j in range(Nw):
            #   create xx,yy,xy   
            y1=np.sum(Nsubap[:j])*2
            y2=y1+Nsubap[j]
            sizeY=Nsubap[j]
            dy=sizeY/subdiv
            if(subdiv*dy<sizeY):
                print "not an integer division"
            idy=j*2*subdiv
            idy2=(j*2+1)*subdiv

            xx=Cmm[x1:x1+sizeX, y1:y1+sizeY]
            yy=Cmm[x2:x2+sizeX, y2:y2+sizeY]
            xy=Cmm[x2:x2+sizeX, y1:y1+sizeY]

            #   decompose xx,yy,xy
            for si in range(subdiv):
                for sj in range(subdiv):
                    blocks[idx+si ,idy+sj ]=xx[si*dx:(si+1)*dx,sj*dy:(sj+1)*dy]
                    blocks[idx2+si,idy2+sj]=yy[si*dx:(si+1)*dx,sj*dy:(sj+1)*dy]
                    blocks[idx+si ,idy2+sj]=xy[si*dx:(si+1)*dx,sj*dy:(sj+1)*dy]
                    blocks[idx2+si,idy+sj ]=xy[si*dx:(si+1)*dx,sj*dy:(sj+1)*dy].T


    return [blocks,size]#,sizex,sizey]


def svd_block(blocks):
    """return an array containing for each position the eigenvalues of
    the corresponding position in the input
    """
    eigenV=np.copy(blocks)
    for i in range(eigenV.shape[0]):
        for j in range(eigenV.shape[1]):
            u,eigenV[i,j],v=np.linalg.svd(eigenV[i,j])

    return eigenV


def  get_nbEigenval(eig,thresh=0,norm=1):
    """return the number of eigenvalues to keep given a threshold
    """
    nbEig=np.copy(eig)
    for i in range(eig.size):
        if(norm==1):
            nbEig.itemset(i,float(np.sum(nbEig.item(i)>thresh))/eig.item(i).size)
        else:
            nbEig.itemset(i,np.sum(nbEig.item(i)>thresh))

    return nbEig.astype(np.float64)


def get_Nsubap(tomo):
    Nsubap=np.zeros((2*tomo.Nw,2*tomo.Nw))
    for i in range(tomo_target.Nw):
        for j in range(tomo_target.Nw):
            Nsubap[2*i,2*j]=min(tomo_target.Nsubap[i],tomo_target.Nsubap[j])
            Nsubap[2*i+1,2*j]=Nsubap[2*i,2*j]
            Nsubap[2*i,2*j+1]=Nsubap[2*i,2*j]
            Nsubap[2*i+1,2*j+1]=Nsubap[2*i,2*j]

    return Nsubap


def plot_nbEig(nbEig,size=None):
    """plot the matrix nbEig (with a colorbar)

        if size is not None: resize the nbEig matrix with blocks of size 'size'
    """
    if(size==None):
        M=nbEig
    else:
        end  =np.cumsum(size)
        start=end-size
        M=np.zeros((end[-1],end[-1]))
        for i in range(size.size):
            for j in range(size.size):
                M[start[i]:end[i], start[j]:end[j]]=nbEig[i,j]

    fig=pl.figure()
    ax=fig.add_subplot(111)
    ax.set_title("percentage of eigenvalues kept")
    cax=ax.matshow(M)
    fig.colorbar(cax)
    pl.show()

def doAll(tomo,cmm=None,subdiv=1,thresh=0,norm=1,scale=1):
    """return an array containing the array of eigenvalues for each blocks 
    """
    data=all_block(tomo,Cmm=cmm,subdiv=subdiv)
    data[0]=svd_block(data[0])
    data_p=get_nbEigenval(data[0],thresh=thresh,norm=norm)
    if(scale==1):
        plot_nbEig(data_p,data[1])
    else:
        plot_nbEig(data_p)
    return data

# draw an other matrix with a different threshold
def tryAgain(data,thresh=0,size=None,norm=1):
    data=get_nbEigenval(data,thresh=thresh,norm=norm)
    plot_nbEig(data,size)



# memory footprint 
def get_size(data,thresh=0):
    storage=0
    N=data[1].size
    for i in range(N):
        x=data[1][i]
        for j in range(N):
            y=data[1][j]

            k=np.sum(data[0][i,j]>thresh)
            if(k<x*y/(x+y)):
                storage+=( x+x)*k
            else:
                storage+=data[1][i]*data[1][j]

    return float(8*storage)

# ratio memory footprint compressed/fully stored 
def print_compRate(data,cmm,list_thresh):
    store_full=float(8*cmm.size)
    for i in list_thresh:
        print "storage rate for a threshold of", i,":",get_size(data,thresh=i)/store_full

