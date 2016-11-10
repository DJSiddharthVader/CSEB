import numpy as np
import numpy.random as npr
import matplotlib.pyplot as plt

#No Vital Dynamics

def infectionrate(S0, probI):
    """takes a probability(probS) in decimal form and finds the distribution of states over TotalPop
    random trials (the chain binomial), return the number of S individuals at a given infection prob
    after infection is allowed to occur"""
    inf = 0
    for x in range(S0):
        bi = npr.uniform(0, 1)
        if bi < probI:
            inf += 1
        else:
            inf += 0
    return inf

def recoveryrate(I0, probR):
    """ using the amount of infected it calculates the amount of people that recover from
    the infection
    :param IR: infected population
    :param probR: recovery probability
    :return: recovered individuals
    """
    rec = 0
    for x in range(I0):
        bi = npr.uniform(0, 1)
        if bi < probR:
            rec += 1
        else:
            rec += 0
    return rec

def cb_gen(S0, I0, R0, probInf, probRe):
    """find the number of infected and succeptible individuals in a population over 1 generation under
    the chain binomial epidemic model given a probability of infection(prob), the initial number of
     succeptible individuals(S)and the initial number of infected individuals (I0). It returns
     the new number of succeptible and infected individuals."""

    PI = 1- (1 - probInf) ** (I0)
    new_S = S0 - infectionrate(S0, PI)
    #new_S = S0 - infectionrate(S0, probInf)
    new_R = recoveryrate(I0, probRe) + R0
    new_I = I0 + S0 - new_S - new_R + R0
    if new_I + new_S + new_R != S0 + I0 + R0:
        raise ValueError
    return new_S, new_I, new_R

#cb_gen(20, 5, 7, 0.2, 0.4)

def cb_sim(S0, I0, R0, probSuc, probRec, nmax = 10):
    """simulates the chain binomial epidemic model (cb_gen) for nmax generations (10 by deafault) and returns the number of
    infected individuals in each generation in a tuple (IG)"""
    SG = []
    IG = []
    RG = []
    SG.append(S0)
    IG.append(I0)
    RG.append(R0)
    Snew, Inew, Rnew = cb_gen(S0, I0, R0, probSuc, probRec)
    SG.append(Snew)
    IG.append(Inew)
    RG.append(Rnew)
    for i in range(nmax - 1):
        Snew, Inew, Rnew = cb_gen(Snew, Inew, Rnew, probSuc, probRec)
        if Inew == 0:
            SG.append(Snew)
            IG.append(Inew)
            RG.append(Rnew)
            return tuple(SG), tuple(IG), tuple(RG)
        SG.append(Snew)
        IG.append(Inew)
        RG.append(Rnew)
    return tuple(SG), tuple(IG), tuple(RG)

def cb_sim_graph(S0 = 2192, I0 = 3, R0 = 5, probInf = 0.01, probRec = 0.02, nmax = 20):
    '''Uses the cb_sim() function to generate data for one instance of an infection (either 10 generations or until there
    are 0 infected people. It then plots the number of succeptibles, infected and recovered against the number of generations.
    It also prints the Y vlaues for each populations
    '''
    x, y, z, = cb_sim(S0, I0, R0, probInf, probRec, nmax)
    print("Susceptible", x,'\n',"Infected  ", y, '\n' "Recovered  ",z)
    plt.plot( x, 'r--', y, 'g-', z, 'b:')
    plt.axis([0,len(x),-1,S0+I0+R0])
    totpop = S0+I0+R0
    plt.title("SIR Model of a Population of Size %i" %totpop)
    plt.legend(("Succeptible", "Infected", "Recovered"))
    plt.xlabel("Generations")
    plt.ylabel("Populations")
    plt.show()

#cb_sim_graph()

def update_cb_dict(d,k):
    """updates a dictionary d. If key k is already a key, adds one to d[k], if not initializes d[k] = 1"""
    #for z in (k):
    if k in d:
        d[(k)] += 1
    else:
        d[(k)] = 1



def run10kM4(succeptible, infected, recovered, recoveryprob, infectionprobmax, infectionprobstep, gens,sims):
    """ Runs 10000 simulations of the binomial epidemic model using cb_sim() and updates
     the dictionaries S1, I1, R1 (using update_cb_dict) with each simulation's tuples of succeptibles, infected and
     recovered as a key and the number of occurences as the value. Then the average number of individals for each class
      of individuals"""
    probiterate = infectionprobstep
    averageepidemicsizes = []
    while probiterate <= infectionprobmax+infectionprobstep:
        succeptiblegenerations = {}
        infectedgenerations    = {}
        recoveredgenerations   = {}
        for n in range(sims):
            a, b, c = cb_sim(succeptible, infected, recovered, probiterate, recoveryprob, gens)
            update_cb_dict(succeptiblegenerations, a)
            update_cb_dict(infectedgenerations, b)
            update_cb_dict(recoveredgenerations, c)

        individualepidemicsize = 0
        # print("infected dict", infectedgenerations)
        # print("recoverd dict", recoveredgenerations)
        for key in recoveredgenerations.keys():
            individualepidemicsize += (key[-1] - key[0] - infected)*recoveredgenerations[key]

        averageepidemicsizes.append(individualepidemicsize/sims)
        probiterate += infectionprobstep

    averageepidemicsizes = np.array(averageepidemicsizes)
    print("array of the average epidemic size for each infection probability", '\n', averageepidemicsizes, '\n', len(averageepidemicsizes))
    return averageepidemicsizes


def run10kM4plot(succeptible, infected, recovered, recoveryprob, infectionprobmax, infectionprobstep, gens,sims):
    infectionprobs = np.arange(infectionprobstep, infectionprobmax+infectionprobstep, infectionprobstep)
    print("Infection probabilities", '\n', infectionprobs, '\n', len(infectionprobs))
    epidemicsizes = run10kM4(succeptible, infected, recovered, recoveryprob, infectionprobmax, infectionprobstep, gens,sims)
    # fig = plt.figure(1)
    # ax = fig.add_subplot(111)
    plt.bar(infectionprobs, epidemicsizes)
    plt.axis([infectionprobstep,infectionprobmax+infectionprobstep, min(epidemicsizes),max(epidemicsizes)])
    plt.title("SIR Average Epidemic Size for Multiple Infection Probabilities")
    #ax.set_xticks(infectionprobs)
    plt.xlabel("Infection Probability")
    plt.ylabel("Average Epidemic Size Over %i Simulations" %sims)
    plt.show()

run10kM4plot(90,5,5,0.2,0.5,0.02, 10,1000)


# def run10kplotM2(Suc = 20, Inf = 4, Rec = 4, probI = 0.3, probR = 0.3, nmax = 10):
#     MeanlistS = []
#     MeanlistI = []
#     MeanlistR = []
#     tempS, tempI, tempR = 0, 0, 0
#     probiterate = 0
#     while probiterate <= 1.00:
#         for x in range(10000):
#             tempresults = [cb_sim(Suc, Inf, Rec, probiterate, probR, nmax)]
#             tempresultsS = tempresults[0][0]
#             tempresultsI = tempresults[0][1]
#             tempresultsR = tempresults[0][2]
#             n = len(tempresults[0][0])
#             for i in range(n):
#                 tempS += tempresultsS[i]
#                 tempI += tempresultsI[i]
#                 tempR += tempresultsR[i]
#         MeanlistS.append(tempS/10000)
#         MeanlistI.append(tempI/10000)
#         MeanlistR.append(tempR/10000)
#         temps, tempI, tempR = 0, 0, 0
#         probiterate += 0.05
#     SSS = np.array(MeanlistS)
#     III = np.array(MeanlistI)
#     RRR = np.array(MeanlistR)
#     print(SSS, '\n', III, '\n', RRR)
#     plt.show()
#     fig.savefig("1MP3Assg6.png")
#     sMeanlistS = []
#     iMeanlistI = []
#     rMeanlistR = []
#     tempS, tempI, tempR = 0, 0, 0
#     probiterate = 0
#     while probiterate <= 1.00:
#         for x in range(10000):
#             tempresults = [cb_sim(Suc, Inf, Rec, probI, probiterate, nmax)]
#             tempresultsS = tempresults[0][0]
#             tempresultsI = tempresults[0][1]
#             tempresultsR = tempresults[0][2]
#             n = len(tempresults[0][0])
#             for i in range(n):
#                 tempS += tempresultsS[i]
#                 tempI += tempresultsI[i]
#                 tempR += tempresultsR[i]
#         sMeanlistS.append(tempS/10000)
#         iMeanlistI.append(tempI/10000)
#         rMeanlistR.append(tempR/10000)
#         temps, tempI, tempR = 0, 0, 0
#         probiterate += 0.05
#     a = np.array(sMeanlistS)
#     b = np.array(iMeanlistI)
#     c = np.array(rMeanlistR)
#     print(a, '\n', b, '\n', c)
#     fig = plt.figure(1)
#     fig2 = plt.figure(2)
#     prob = np.arange(0, 1, 0.05)
#     print(len(x), len(y), len(z), len(a), len(b), len(c))
#     ax1 = fig.add_subplot(1, 1, 1)
#     plt.plot
#     g = [1,2,3,4,5]
#     G = np.array(g)
#     plt.plot(G, G)
#     plt.plot(prob,SSS)
#     plt.plot(prob, SSS, 'r--', prob, III, 'bo', prob, RRR, 'g^')
#     ax1.plt.xlabel("Infection Probability")
#     ax2 = fig2.add_subplot(1, 1, 1)
#     plt.plot(prob, a, 'r--', prob, b, 'bo', prob, c, 'g^')
#     ax2.set_xlabel("Recovery Probability")
#     plt.show()

#
# def run10kM3(Suc = 40, Inf = 2, Rec = 4, ProbI = 0.1, ProbR =0.2, nmax = 10):
#     """ Runs 10000 simulations of the binomial epidemic model using cb_sim() and updates
#      the dictionaries S1, I1, R1 (using update_cb_dict) with each simulation's tuples of succeptibles, infected and
#      recovered as a key and the number of occurences as the value. Then the average number of individals for each class
#       of individuals"""
#     probiterate = 0
#     MeanlistS = []
#     MeanlistI = []
#     MeanlistR = []
#
#     while probiterate <= 1.0:
#         S1 = {}
#         I1 = {}
#         R1 = {}
#         lens, leni, lenr = 0, 0, 0
#         for x in range(10000):
#             a, b, c = cb_sim(Suc, Inf, Rec, probiterate,ProbR, nmax)
#             update_cb_dict(S1, a)
#             update_cb_dict(I1, b)
#             update_cb_dict(R1, c)
#
#         temp = 0
#         lenses = [lens, leni, lenr]
#         dicts = [S1, I1, R1]
#         meanlists = [MeanlistS, MeanlistI, MeanlistR]
#         for t in range(3):
#             u = dicts[t]
#             for x in u.keys():
#                 for y in x:
#                     temp += y*u[x]
#                     lenses[t] += u[x]
#             meanlists[t].append(temp / lenses[t])
#
#         SS = np.array(MeanlistS)
#         II = np.array(MeanlistI)
#         RR = np.array(MeanlistR)
#         probiterate += 0.05
#
#     print(SS, II, RR)
#     return SS, II, RR
