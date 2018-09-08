import random
import operator

mutationRate = 0.1
fitnessReq = 0.1 # from 0 to 1 of how the tolerancy of what percentile is acceptable
neuron = [4,4,2,1] # first number in this list should always match the number of inputs, and the last the number of outputs
population = []
evaluation = []
populationSize = 100
generationNum = 10000


InOut = [
    [[0,0,0,0],[0]],
    [[0,0,0,1],[0]],
    [[0,0,1,0],[0]],
    [[0,1,0,0],[0]],
    [[1,0,0,0],[0]],
    [[0,0,1,1],[0]],
    [[0,1,0,1],[0]],
    [[1,0,0,1],[1]],
    [[0,1,1,0],[1]],
    [[1,0,1,0],[0]],
    [[1,1,0,0],[0]],
    [[1,1,1,0],[0]],
    [[1,1,0,1],[0]],
    [[1,0,1,1],[0]],
    [[0,1,1,1],[0]],
    [[1,1,1,1],[0]], #all possible combinations
    
    [[1,0,0,1],[1]], # repeating the number of combinations that return [1] until the ratio of [1] returning and [0] returning are the same
    [[0,1,1,0],[1]],
    [[1,0,0,1],[1]],
    [[0,1,1,0],[1]],
    [[1,0,0,1],[1]],
    [[0,1,1,0],[1]],
    [[1,0,0,1],[1]],
    [[0,1,1,0],[1]],
    [[1,0,0,1],[1]],
    [[0,1,1,0],[1]],
    [[1,0,0,1],[1]],
    [[0,1,1,0],[1]],
    ]

def createNewGene():
    global neuron
    global genes
    temp = []
    
    for layer in range(0,len(neuron)-1):
        temp.append([])
        for neurR1 in range(0,neuron[layer]):
            for neurR2 in range(0,neuron[layer+1]):
                temp[layer].append(float('%.2f'%(random.uniform(-1,1))))
    Bias_Neurons = []
    for biasNeurLayer in range(0,len(neuron)-1):
        Bias_Neurons.append([])
        for biasNeur in range(0,neuron[biasNeurLayer+1]):
            Bias_Neurons[biasNeurLayer].append(float('%.2f'%(random.uniform(-1,1))))
    temp += Bias_Neurons
    return temp


def createGenes(Strain): # aka create mutated baby
    prevStrain = Strain
    for layer in range(0,len(prevStrain)):
        for gen in range(0,len(prevStrain[layer])):
            test = random.random() < mutationRate
            if test:
                prevStrain[layer][gen] = float('%.2f'%(random.uniform(-1,1)))
    return prevStrain

def ICTBI(inpt): #Int convert to Bool Int
    if(inpt > 0):
        return 1
    return 0
        

def decision(Strain,inpt):
    prevLayer = inpt
    curLayer = []
    for layer in range(0,len(neuron)-1): # Every layer
        for SetCurLayerSize in range(0,neuron[layer+1]):
            curLayer.append(0)
        for gene in range(0,len(Strain[layer])):
            preNeur = gene // neuron[layer]
            if prevLayer[preNeur] == 1: #checks to see if it should bother calculating that neuron
                postNeur = gene % neuron[layer+1]
                curLayer[postNeur] += Strain[layer][gene]*prevLayer[preNeur]
        for bias in range(0,len(curLayer)):
            curLayer[bias] += Strain[layer+len(neuron)-1][bias]
            
        for val in range(0,len(curLayer)):
            curLayer[val] = ICTBI(curLayer[val])
            
        prevLayer = curLayer
        curLayer = []
    return prevLayer
                
        
def main():
    for ins in range(0,populationSize):    
        population.append(createNewGene()) #make the babies
    
    for num in range(0,generationNum):
        generation()
        print(evaluation)
    
    
def generation():
    global population
    global evaluation
    for ins in range(0,populationSize): #fill in population gaps with mutated instances
        if population[ins] == 0:
            check = ins
            parentNotFound = True
            while parentNotFound:
                
                if population[check] == 0:
                    check = random.randint(0,populationSize-1)
                else:
                    population[ins] = createGenes(population[check])
                    parentNotFound = False
          
    evaluation = []
    for ins in range(0,populationSize): # evaluate each member in the population to see how effective it is
        evaluation.append([fittness(ins,InOut),ins])
        
    evaluation = sorted(evaluation, key=operator.itemgetter(0),reverse=True)
    
    
    
    for cull in range(round( populationSize * (fitnessReq)),populationSize):
        population[evaluation[cull][1]] = 0
            
            
def fittness(ins,InOut):
    value = 0
    for Q in range(0,len(InOut)):
        output = decision(population[ins],InOut[Q][0])
        if output == InOut[Q][1]:
            value += 1
    return value / len(InOut)

main()

