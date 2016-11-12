import numpy as np
import numpy.random as npr
import matplotlib.pyplot as plt

#No Vital Dynamics

def infectionrate(S0, probI):
    ''' using the amount of infected and the recovery probability it calculates the amount of people that recover from
    the infection
    :param S0: succeptible population
    :param probI: infection probability
    :return: integer of infected individuals
    '''
    if probI < 0:
        raise ValueError("Cannot have negative probability")
    inf = 0
    for x in range(S0):
        bi = npr.uniform(0, 1)
        if bi < probI:
            inf += 1
        else:
            inf += 0
    return inf

def recoveryrate(I0, probR):
    ''' using the amount of infected and the recovery probability it calculates the amount of people that recover from
    the infection
    :param I0: infected population
    :param probR: recovery probability
    :return: integer of recovered individuals
    '''
    if probR < 0:
        raise ValueError("Cannot have negative probability")
    rec = 0
    for x in range(I0):
        bi = npr.uniform(0, 1)
        if bi < probR:
            rec += 1
        else:
            rec += 0
    return rec











def cb_gen_nochain(S0, I0, R0, probInf, probRe):
    '''find the number of infected and succeptible individuals in a population over 1 generation under
    the chain binomial epidemic model given a probability of infection(prob), the initial number of
     succeptible individuals(S)and the initial number of infected individuals (I0). It returns
     the new number of succeptible and infected individuals.
    :param S0:  Initial number of succeptibles individuals
    :param I0: Initial number of infected individuals
    :param R0: Initial number of recovered individuals
    :param probInf: probability that someone will move from the succeptible class to the infected class
    :param probRe: probability that someone will move from the infected class to the recovered class
    :return: three integers, the number of succeptibles, infected and recovered respectively
    '''
    new_S = S0 - infectionrate(S0, probInf)
    new_R = recoveryrate(I0, probRe) + R0
    new_I = I0 + S0 - new_S - new_R + R0
    if new_I + new_S + new_R != S0 + I0 + R0:
        raise ValueError
    return new_S, new_I, new_R

def cb_sim_nochain(S0, I0, R0, probInf, probRec, nmax):
    '''runs a simulation of nmax generation with constant infection probability using cb_gen_nochain()
    :param S0: Initial number of succeptibles individuals
    :param I0: Initial number of infected individuals
    :param R0: Initial number of recovered individuals
    :param probInf: probability that someone will move from the succeptible class to the infected class
    :param probRec: probability that someone will move from the infected class to the recovered class
    :param nmax: the number of generations to run, if the number of infected does not reach zero beforehand
    :return: three tuples, each containing the number of individuals in each class at each generation
    '''
    SG = []
    IG = []
    RG = []
    SG.append(S0)
    IG.append(I0)
    RG.append(R0)
    Snew, Inew, Rnew = cb_gen_nochain(S0, I0, R0, probInf, probRec)
    SG.append(Snew)
    IG.append(Inew)
    RG.append(Rnew)
    for i in range(nmax - 1):
        Snew, Inew, Rnew = cb_gen_nochain(Snew, Inew, Rnew, probInf, probRec)
        if Inew == 0:
            SG.append(Snew)
            IG.append(Inew)
            RG.append(Rnew)
            return tuple(SG), tuple(IG), tuple(RG)
        SG.append(Snew)
        IG.append(Inew)
        RG.append(Rnew)
    return tuple(SG), tuple(IG), tuple(RG)

def cb_sim_graph_nochain(S0 = 1997, I0 = 3, R0 = 0, probInf = 0.01, probRec = 0.02, nmax = 400):
    '''Uses the cb_sim_nochain() function to generate data for one instance of an infection (either 400 generations or until there
    are 0 infected people. It then plots the percentage of succeptibles, infected and recovered against the number of generations.
    :param S0: Initial number of succeptibles individuals
    :param I0: Initial number of infected individuals
    :param R0: Initial number of recovered individuals
    :param probInf: probability that someone will move from the succeptible class to the infected class
    :param probRec: probability that someone will move from the infected class to the recovered class
    :param nmax: the number of generations to run, if the number of infected does not reach zero beforehand
    :return: a graph of the percentage of the population in each class over time
    '''
    susc, infs, recov, = cb_sim_nochain(S0, I0, R0, probInf, probRec, nmax)
    resultlist = [susc, infs, recov]
    suscnew, infsnew, recovnew = [], [], []
    news = [suscnew, infsnew, recovnew]
    for new in range(3):
        resultselect = resultlist[new]
        newsselected = news[new]
        for element in resultselect:
            newsselected.append((element/(S0+I0+R0))*100)
    if len(suscnew) != len(infsnew) or len(suscnew) != len(recovnew) or len(infsnew) != len(recovnew):
        raise ValueError("Classes are of different generation sizes")
    plt.plot(suscnew,'r--', infsnew, 'g-', recovnew, 'b:')
    plt.axis([0,len(suscnew),0,100])
    totpop = S0+I0+R0
    plt.title("SIR Model of a Population of Size %i" %totpop)
    plt.legend(("Succeptible", "Infected", "Recovered"))
    plt.xlabel("Generations")
    plt.ylabel("Percentage of Population in Each Class")
    plt.show()

cb_sim_graph_nochain()







def cb_gen(S0, I0, R0, probInf, probRe):
    '''find the number of infected and succeptible individuals in a population over 1 generation under
    the chain binomial epidemic model, where infection probability vaires with the number of infected individuals
    based on the equation:
    Infection Probability = 1 - (1 - Infection Probability)^(initial number of infected)
    :param S0:  Initial number of succeptibles individuals
    :param I0: Initial number of infected individuals
    :param R0: Initial number of recovered individuals
    :param probInf: probability that someone will move from the succeptible class to the infected class
    :param probRe: probability that someone will move from the infected class to the recovered class
    :return: three integers, the number of succeptibles, infected and recovered respectively
    '''
    PI = 1- (1 - probInf) ** (I0)
    new_S = S0 - infectionrate(S0, PI)
    new_R = recoveryrate(I0, probRe) + R0
    new_I = I0 + S0 - new_S - new_R + R0
    if new_I + new_S + new_R != S0 + I0 + R0:
        raise ValueError
    return new_S, new_I, new_R

def cb_sim(S0, I0, R0, probSuc, probRec, nmax):
    '''runs a simulation of nmax generation with constant infection probability using cb_gen()
    :param S0: Initial number of succeptibles individuals
    :param I0: Initial number of infected individuals
    :param R0: Initial number of recovered individuals
    :param probInf: probability that someone will move from the succeptible class to the infected class
    :param probRec: probability that someone will move from the infected class to the recovered class
    :param nmax: the number of generations to run, if the number of infected does not reach zero beforehand
    :return: three tuples, each containing the number of individuals in each class at each generation
    '''
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

def cb_sim_graph(S0 = 1997, I0 = 3, R0 = 0, probInf = 0.01, probRec = 0.02, nmax = 20):
    '''Uses the cb_sim() function to generate data for one instance of an infection (either 400 generations
    or until there are 0 infected people. It then plots the percentage of succeptibles, infected and recovered
    against the number of generations.
    :param S0: Initial number of succeptibles individuals
    :param I0: Initial number of infected individuals
    :param R0: Initial number of recovered individuals
    :param probInf: probability that someone will move from the succeptible class to the infected class
    :param probRec: probability that someone will move from the infected class to the recovered class
    :param nmax: the number of generations to run, if the number of infected does not reach zero beforehand
    :return: a graph of the percentage of the population in each class over time
    '''
    susc, infs, recov, = cb_sim(S0, I0, R0, probInf, probRec, nmax)
    resultlist = [susc, infs, recov]
    suscnew, infsnew, recovnew = [], [], []
    news = [suscnew, infsnew, recovnew]
    for new in range(3):
        resultselect = resultlist[new]
        newsselected = news[new]
        for element in resultselect:
            newsselected.append((element/(S0+I0+R0))*100)
    if len(suscnew) != len(infsnew) or len(suscnew) != len(recovnew) or len(infsnew) != len(recovnew):
        raise ValueError("Classes are of different generation sizes")
    plt.plot(suscnew,'r--', infsnew, 'g-', recovnew, 'b:')
    plt.axis([0,len(suscnew),0,100])
    totpop = S0+I0+R0
    plt.title("SIR Model of a Population of Size %i" %totpop)
    plt.legend(("Succeptible", "Infected", "Recovered"))
    plt.xlabel("Generations")
    plt.ylabel("Percentage of Population in Each Class")
    plt.show()

cb_sim_graph()










def update_cb_dict(d,k):
    '''updates a dictionary d. If key k is already a key, adds one to d[k], if not initializes d[k] = 1
    :param d: dictionary
    :param k: iterable
    :return: dictionary with keys that are elements of k
    '''

    if k in d:
        d[(k)] += 1
    else:
        d[(k)] = 1

def run10kM4(succeptible, infected, recovered, recoveryprob, infectionprobmax, infectionprobstep, gens,sims):
    '''
    :param succeptible: initial number of succeptible individuals
    :param infected: initial number of infected individuals
    :param recovered: initial number of recovered individuals
    :param recoveryprob: probability of recovery for an infected individuals
    :param infectionprobmax: the maximum infection probability being tested
    :param infectionprobstep: how much the infection probabilty is being incremented
    :param gens: maximum number of generations being run
    :param sims: number of simulations to run at each infection probability
    :return: array of the average epidemic size over sims simulations at each infection proability
    '''
    probiterate = infectionprobstep
    averageepidemicsizes = []
    while probiterate <= infectionprobmax:
        succeptiblegenerations = {}
        infectedgenerations  = {}
        recoveredgenerations = {}
        for n in range(sims):
            a, b, c = cb_sim(succeptible, infected, recovered, probiterate, recoveryprob, gens)
            update_cb_dict(succeptiblegenerations, a)
            update_cb_dict(infectedgenerations, b)
            update_cb_dict(recoveredgenerations, c)
        # print(succeptiblegenerations, '\n', max(succeptiblegenerations.values()))

        individualepidemicsize = 0
        for key in recoveredgenerations.keys():
            individualepidemicsize += (key[-1] - key[0]-infected)*(recoveredgenerations[key]/sims)
        averageepidemicsizes.append(individualepidemicsize)

        # individualepidemicsize , recoveredsize = 0, 0
        # for key in recoveredgenerations.keys():
        #     recoveredsize += key[0]
        #     recoveredsize -= key[-1]
        #     recoveredsize -= infected
        #     recoveredsize /= succeptible
        #     recoveredsize *= recoveredgenerations[key]
        #     recoveredsize /= sims
        #     individualepidemicsize += recoveredsize
        #     recoveredsize = 0
        #
        # averageepidemicsizes.append((individualepidemicsize))
        probiterate += infectionprobstep

    averageepidemicsizes = np.array(averageepidemicsizes)
    print("array of the average epidemic size for each infection probability", '\n', averageepidemicsizes, '\n', len(averageepidemicsizes))
    return averageepidemicsizes


def run10kM4plot(succeptible, infected, recovered, recoveryprob, infectionprobmax, infectionprobstep, gens,sims):
    '''
    :param succeptible: initial number of succeptible individuals
    :param infected: initial number of infected individuals
    :param recovered: initial number of recovered individuals
    :param recoveryprob: probability of recovery for an infected individuals
    :param infectionprobmax: the maximum infection probability being tested
    :param infectionprobstep: how much the infection probabilty is being incremented
    :param gens: maximum number of generations being run
    :param sims: number of simulations to run at each infection probability
    :return: scatter plot and bar graph of the average epidemic size for each infection probability
    '''
    if (infectionprobmax*10e5)%(infectionprobstep*10e5) != 0:
        print(infectionprobmax%infectionprobstep)
        raise ValueError("cannot step through probabilites properly, make sure infectionprobstep equally divides infectionprobmax")

    infectionprobs = np.arange(infectionprobstep, infectionprobmax+infectionprobstep/2, infectionprobstep)
    epidemicsizes = run10kM4(succeptible, infected, recovered, recoveryprob, infectionprobmax+infectionprobstep/2, infectionprobstep, gens,sims)
    print("Infection probabilities", '\n', infectionprobs, '\n', len(infectionprobs))
    tenpercentgraphspace = (max(epidemicsizes) - min(epidemicsizes))/10

    plt.figure(1)
    plt.scatter(infectionprobs, epidemicsizes)
    plt.axis([0,infectionprobmax+infectionprobstep, min(epidemicsizes) - tenpercentgraphspace,max(epidemicsizes)+tenpercentgraphspace])
    plt.xticks(infectionprobs)
    plt.title("SIR Average Epidemic Size for Multiple Infection Probabilities")
    plt.xlabel("Infection Probability")
    plt.ylabel("Average Epidemic Size Over %i Simulations" %sims)

    plt.figure(2)
    plt.bar(infectionprobs, epidemicsizes, width=infectionprobstep)
    plt.axis([infectionprobstep,infectionprobmax+infectionprobstep, min(epidemicsizes) - tenpercentgraphspace,max(epidemicsizes)+tenpercentgraphspace])
    plt.xticks(infectionprobs)
    plt.title("SIR Average Epidemic Size for Multiple Infection Probabilities")
    plt.xlabel("Infection Probability")
    plt.ylabel("Average Epidemic Size Over %i Simulations" %sims)
    plt.show()

run10kM4plot(400,5,5,0.3,0.5,0.02, 10,100)


