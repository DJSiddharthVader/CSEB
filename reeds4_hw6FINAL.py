import numpy as np
import numpy.random as npr
import matplotlib.pyplot as plt

def rbinom(N,prob):
    """takes a probability(prob) in decimal form and finds the distribution of states over N
    random trials (the chain binomial)"""
    Phy = 0
    for x in range(N):
        bi = npr.uniform(0, 1)
        states = []
        if bi > prob:
            states.extend('clean')
        else:
            states.append('infected')
        for x in states:
            if x == "infected":
                Phy += 1
    return Phy



def cb_gen(S, I0, prob):
    """find the number of infected and succeptible individuals in a population over 1 generation under
    the chain binomial epidemic model given a probability of infection(prob), the initial number of
     succeptible individuals(S)and the initial number of infected individuals (I0). It returns
     the new number of succeptible and infected individuals."""
    P = 1-(1-prob)**(I0)
    new_I = rbinom(S, P)
    new_S = S - new_I
    return (new_S, new_I)




def cb_sim(S, I0, prob, nmax = 10):
    """simulates the chain binomial epidemic model for nmax generations (10 by deafault) and returns the number of
    infected individuals in each generation in a tuple (IG)"""
    IG = []
    IG.append(I0)
    Snew, Inew = cb_gen(S, I0, prob)
    IG.append(Inew)
    for i in range(nmax):
        Snew, Inew = cb_gen(Snew, Inew, prob)
        if Inew == 0:
            return tuple(IG)
        IG.append(Inew)
    return tuple(IG)




def update_cb_dict(d,k):
    """updates a dictionary d. If key k is already a key, adds one to d[k], if not initializes d[k]"""
    #for z in (k):
    if k in d:
        d[(k)] += 1
    else:
        d[(k)] = 1


        
        
def run10k():
    """ Runs 10000 simulations of the binomial epidemic model using cb_sim(5,2,0.4,10) and updates
     the dictionary EO with each simulation's infection as a key and the number of occurences as the
     value. Then prints the most likely infection outcome (key with largest value).
    """
    E1 = {}
    for x in range(10000):
        y = cb_sim(5, 1, 0.4, 10)
        update_cb_dict(E1, y)
    print (E1)
    for z in E1.keys():
        if E1[z] == max(E1.values()):
            print("Most likely outcome:", z)
        else:
            pass
    return (E1)




def run10kplotM2():
    '''Runs 10000 simulations of the binomial epidemic model using cb_sim(5,2,0.4,10) for multiple infection rates (0 to 1 with a 0.05 
    step value). creates a grap with infection probability on the x axis and the mean number of infected for the simulation of 
    each probability'''
    Meanlist = []
    temp = 0
    z = 0
    probs = np.arange(0, 1, 0.05)
    while z <= 1.00:
        for x in range(10000):
            for x in cb_sim(5, 1, z, 10):
                temp += x
        Meanlist.append(temp/10000)
        temp = 0
        z += 0.05
    y = np.array(Meanlist)
    print(y)
    fig = plt.figure()
    ax1 = fig.add_subplot(1, 1, 1)
    plt.scatter(probs, y)
    ax1.set_xlabel("Infection Probability")
    ax1.set_ylabel("Mean Infection Number")
    fig.suptitle("Mean Infection Number vs. Infection Probability")
    plt.show()
    fig.savefig("1MP3Assg6.png")

