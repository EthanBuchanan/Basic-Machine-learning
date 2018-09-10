import random
import operator
import copy
import xml.etree.ElementTree as ET
from xml.etree.ElementTree import Element, SubElement, Comment, ElementTree

# User changing variables (change as needed)


fitnessReq = 0.1 # from 0 to 1 of how the tolerancy of what percentile is acceptable
mutationRate = 0.5
populationSize = 100
generationNum = 500
neuron = [4,4,2,1] # first number in this list should always match the number of inputs, and the last the number of outputs



file_Path = "XML_File_Test.xml" # the path to the file where generations are stored
#make sure that if you change the values in the 'neuron' list that you also change the file_Path 

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

# declaration variables (don't touch)

population = []
evaluation = []
try:
    tree = ET.parse(file_Path)
except:
    test = """
<data>
</data>
"""
    tree = ElementTree(ET.fromstring(test))
    tree.write(file_Path)

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

def loadPrev():
    global tree
    root = tree.getroot()
    
    instances = populationSize
    if len(root) < populationSize:
        instances = len(root)
    for saved in range(instances):
        population[saved] = []
        for layer in range(len(root[saved])):
            entry = root[saved][layer].text
            entry = entry.split()
            entryL = list(map(float,entry))
            population[saved].append(entryL)
    
def saveCur():
    test = """
<data>
</data>
"""
    tree = ElementTree(ET.fromstring(test))
    root = tree.getroot()
    for saved in range(populationSize):
        if population[saved] != 0:
            ins = ET.Element('ins')
            for layer in range(len(population[saved])):
                elem = ET.SubElement(ins,'layer')
                elem.text = (str(population[saved][layer])[1:-1]).replace(",","")
            root.append(ins)
    
    
    tree.write(file_Path)



    
        
def main():
    for ins in range(0,populationSize):    
        population.append(0) #make the babies
    population[0] = createNewGene()    
    loadPrev()

    for num in range(0,generationNum):
        generation()
        print(num)
    print(evaluation)
    saveCur()

def DuplicateSeparateList(startList):
    newList = []
    for entry in startList:
        newList.append(entry)
    return newList

def generation():
    global population
    global evaluation
    
    for ins in range(populationSize): #fill in population gaps with mutated instances
        if population[ins] == 0:
            check = ins
            parentNotFound = True
            while parentNotFound:
                if population[check] == 0:
                    check = random.randint(0,populationSize-1)
                else:
                    
                    population[ins] = createGenes(copy.deepcopy(population[check]))
                    parentNotFound = False
    
    
    evaluation = []
    for ins in range(0,populationSize): # evaluate each member in the population to see how effective it is
        evaluation.append([fittness(ins),ins])
        
    evaluation = sorted(evaluation, key=operator.itemgetter(0),reverse=True)
    
    
    standard = int( populationSize * (fitnessReq))
    if standard < 1:
        standard = 1
    for cull in range(standard,populationSize):
        population[evaluation[cull][1]] = 0
        evaluation[cull] = 0
    
            
def fittness(ins):
    global InOut
    value = 0
    for Q in range(0,len(InOut)):
        output = decision(population[ins],InOut[Q][0])
        if output == InOut[Q][1]:
            value += 1
    return value / len(InOut)

main()


