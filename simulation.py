from matplotlib import pyplot as plt
import numpy as np
import math
from scipy.linalg import expm




def Lanczos(A,n,b_i):
    m = len(A)

    #Lanczos Iteration
    Q = []
    Beta =[]
    q0 = []
    q1 = []
    for i in range(m):
        q0.append(0)
        q1.append(0)
    q1[b_i]=1
    Q.append(np.array(q0))
    Q.append(np.array(q1))
    Beta.append(0)

    for i in range(n):
        v = np.dot(A, Q[i+1])
        alpha = np.dot(np.transpose(Q[i+1]),v)
        v = np.subtract(v,np.subtract(np.dot(Beta[i],Q[i]),np.dot(alpha,Q[i+1])))
        beta_n = np.linalg.norm(v)
        qn = np.dot(v, (1 / beta_n))
        Beta.append(beta_n)
        Q.append(np.array(qn))
    del Q[0]
    del Q[1]
    Q = np.transpose(Q)
    return np.dot(np.transpose(Q),np.dot(A,Q)), Q



def compute_centrality(n):
    A = np.loadtxt(open("A.csv", "rb"), delimiter=",")

    W = A

    for i in range(len(W)):
        W[i][i]=0
    beta = 0.01
    c= []
    #Define e1
    for j in range(len(A)):
        T = Lanczos(W, n, j)[0]
        e = []
        for i in range(n):
            e.append(0)
        e[0] = 1
        e1 = []
        e1.append(e)
        e1 = np.transpose(e1)
        ci = np.dot(np.transpose(e1),np.dot(expm(np.dot(beta,T)),e1))
        c.append(ci)
    return c

def compute_u(t,n):

    A = np.loadtxt(open("A.csv", "rb"), delimiter=",")

    L = A

    T,Q = Lanczos(L, n, 1)
    e = []
    for i in range(n):
        e.append(0)
    e[0] = 1
    e1 = []
    e1.append(e)
    e1 = np.transpose(e1)

    ui = np.dot(Q,np.dot(expm(np.dot(-t,T)),e1))
    return ui


def compute_p(t,n,l,d):

    I = []
    p = []

    I.append(1)
    u_0 = compute_u(0,n)
    p_0 = []
    for i in u_0:
        p_0.append(abs(min(1,i[0])))
    p.append(p_0)

    for i in range(l):
        In = 0
        p_n = []
        for j in p[i]:
            if j >= d:
                In = In + 1
        u_n = compute_u(t/l,n)
        for k in u_n:
            p_n.append(abs(min(1, In*k)))
        p.append(p_n)
        I.append(In)

    return p, I


def disease_simulation():
    # n = 10
    # c = compute_u(n,1)
    # c_10 = compute_u(n+10,1)
    #
    # while(abs(c-c_10)>=0.0001):
    #     n = n+10
    #     c = compute_u(n,1)
    #     c_10 = compute_u(n + 10,1)
    #     print(n)
    #     print(abs(c - c_10))
    # return n

    c = compute_centrality(70)
    #c = compute_u(30,1)
    #p,I = compute_p(0.2, 30, 100, 0.001)

    #print(I[-1])
    # Y= []
    # for i in range(len(p[-1])):
    #     Y.append(p[-1][i])
    # X = list(range(1, len(p[-1])+1))
    # plt.figure()
    # plt.ylabel('p_i')
    # plt.xlabel('i')
    # plt.plot(X, Y)
    # plt.title("Plotting p_i, where  t = 0.2 and d = 0.001")
    # plt.show()

    Y = []
    X = []
    # for i in range(15):
    #     Y.append(I[i])
    #     X.append((i+1)*(1/100)*(0.2))
    # plt.figure()
    # plt.ylabel('I(t)')
    # plt.xlabel('t')
    # plt.plot(X, Y)
    # plt.title("Plotting I(t)")
    # plt.show()

    for i in range(len(c)):
        Y.append(c[0][0][i])
    X = list(range(1, len(c)+1))
    plt.figure()
    plt.ylabel('c_i')
    plt.xlabel('i')
    plt.plot(X, Y)
    plt.yscale('log')
    plt.title("Plotting c_i, where  β = 0.01")
    plt.show()

    return 1


disease_simulation()



