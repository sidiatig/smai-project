import numpy as np
import sys
from numpy import linalg as LA

def jump_detection(signal):
    signal = signal.T
    L = len(signal)
    M = np.floor(np.log2(L))

    if (M >= 7):
        niter = 4
    elif (M >= 4):
        niter = 2
    else:
        niter = 0

    if (niter == 0):
        ligne_cour = signal
        res = abs(np.diff(ligne_cour, axis=0))
        H = np.zeros((1,L-2))
        temp = res.flatten('F')
        H = temp[0:-1] + temp[1:]
        Ind = []
        for i in range(3, L-3):
            if (temp[i] > max(temp[i+1:i+3])):
                Ind.append(i)
        Indn = Ind
    else:
        nbpoint = 4
        ligne_cour = signal
        label = np.zeros((niter, L))
        for k in range(0, niter):
            L = len(ligne_cour)
            H = np.zeros((1, L-2))
            res = abs(np.diff(ligne_cour, axis=0))
            temp = res.flatten('F')
            H = temp[0:-1] + temp[1:]
            Ind = []
            for i in range(2+3*(nbpoint-k+1)-1, L-2-3*(nbpoint-k+1)):
                if (H[i] > max(H[i+1:i+3*(nbpoint-k+1)+1]) and (H[i-1] > max(H[i-1-3*(nbpoint-k+1):i-1]))):
                    Ind.append(i)
            q = 0
            Indn = []
            while (q < len(Ind)-1):
                if (Ind[q] + 1 == Ind[q+1]):
                    if (H[Ind[q+1]] >= H[Ind[q]]):
                        Indn.append(Ind[q+1])
                    else:
                        Indn.append(Ind[q])
                    q = q+2
                else:
                    Indn.append(Ind[q])
                    q = q+1
            # possible source of errors
            if (q+1 == len(Ind)):
                Indn.append(Ind[q])

            label[k][Indn] = np.ones((1,len(Indn)))
            temp2 = ligne_cour.flatten('F')

            ligne_cour = (temp2[0:-1:2] + temp2[1::2])/2

        # check here if errors found
        Xint = np.arange(1, len(signal)+1)*(label[niter-1])
        Xint = (2**(niter-1))*Xint[Xint>0]
        Xsup = np.arange(1, len(signal)+1)*(label[0])
        Xsup = Xsup[Xsup>0]

        N = len(Xint)
        Ns =len(Xsup)
        Indn = []
        for i in range(0, N):
            for j in range(0, Ns):
                if (abs(Xint[i]-Xsup[j])<2*niter):
                    Indn.append(Xsup[j])

    return Indn

def sinkhorn(M,mu,nu,lam,iterations):
    d=len(M)
    a=np.ones(d)
    tolerance=1e-10
    Err=np.empty([iterations,2])
    k=np.exp(np.multiply(M,-1*lam))
    # print k
    for i in range(0,iterations):
        new_k=np.matmul(np.transpose(k),a)
        b=np.divide(nu,new_k)
        Err[i][0]=LA.norm(a*np.matmul(k,b)-mu,1)
        a=np.divide(mu,np.matmul(k,b))
        Err[i][1]=LA.norm(b*np.matmul(np.transpose(k),a)-nu,1)
        if max(Err[i][0],Err[i][1]) < tolerance:
            break

    gamma=np.matmul(np.diag(a),k)
    gamma=np.matmul(gamma,np.diag(b))
    return [a,b,gamma]

def cluster(idx_edges,idxb):
    idx_edges=np.array(idx_edges).astype(int)
    d=len(idxb)
    l=len(idx_edges)
    clu=np.empty([d+1,1])
    ans=np.empty([d+1,1])
    co=np.zeros([l+2,1])
    clu[0:idx_edges[0]+1]=1
    co[1]=idx_edges[0]+1
    for i in range(1,len(idx_edges)):
        clu[idx_edges[i-1]+1:idx_edges[i]+1]=i+1
        co[i+1]=idx_edges[i]-idx_edges[i-1]
    clu[idx_edges[l-1]+1:d+1]=l+1
    co[l+1]=d-idx_edges[l-1]
    for i in range(0,len(idxb)):
        ans[i]=clu[idxb[i]]

    return [co,ans]

if __name__ == "__main__":
    f=open("./ml-100k/u.item","r")
    movies=[]
    for line in f:
        line=line.split('\n')[0]
        line=line.split('|')
        movies.append(line[1])

    filename = sys.argv[1]
    filereader = open(filename,"r")
    user = []
    movie = []
    rating = []
    for line in filereader:
        line = line.split('\n')[0]
        line = line.split('\t')
        user.append(line[0])
        movie.append(line[1])
        rating.append(line[2])
    nm = 0
    nu = int(user[len(user)-1])
    ns=int(len(user))
    for i in movie:
        if nm <= int(i):
            nm = int(i)
    data = np.zeros((nu,nm))
    for i in range(len(user)):
        u = int(user[i])
        j = int(movie[i])
        data[u-1][j-1] = int(rating[i])



    data_trans=np.transpose(data)
    z=data_trans[np.random.choice(data_trans.shape[0],len(data_trans[0]), replace=True)]
    z=np.transpose(z)
    z_trans=np.transpose(z)
    M=np.empty([z.shape[0],z.shape[1]])
    for i in range(z.shape[0]):
         for j in range(z.shape[1]):
             M[i][j]=np.sum(np.square(z[i]-z_trans[j]))

    cmean=np.count_nonzero(z,axis=1)
    cmean=cmean.astype(float)/ns
    rmean=np.count_nonzero(z,axis=0)
    rmean=rmean.astype(float)/ns

    # print M.max()
    al,be,gam=sinkhorn(M,cmean,rmean,0.001,10)
    idx_al=np.argsort(al)
    idx_be=np.argsort(be)
    L_al = jump_detection(np.sort(al))
    L_be = jump_detection(np.sort(be))

    g=len(L_al)+1
    m=len(L_be)+1

    [r_co,row_cluster_id]=cluster(L_al,idx_al)
    [c_co,col_cluster_id]=cluster(L_be,idx_be)
    avg_movie_rating=np.mean(data,axis=0,dtype=float)
    clusters_num=np.amax(col_cluster_id)
    clusters_num=int(clusters_num)
    max_five=np.empty([int(clusters_num),5])
    a=np.zeros([clusters_num,int(len(data[0]))])
    c=0
    for i in range(clusters_num):
        c=0
        for j in range(col_cluster_id.shape[0]):
            if int(col_cluster_id[j][0])==i+1:
                a[i][c]=j+1
                c=c+1
        np.sort(a[i])[::-1]
    top_movies=[]
    for i in range(clusters_num):
        mov=[]
        for j in range(5):
            mov.append(movies[int(a[i][j])])
        top_movies.append(mov)
    print top_movies
