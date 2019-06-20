import numpy as np


eps = 10e-8

### class definitions
class cell:
    def __init__(self, x, y):
        self.x = x ## x coordinate of grid cell
        self.y = y ## y coordinate of grid cell
        self.contents = None ## pointer to agent in cell
        self.concentration = 0.0 ## the concentration of antibiotic resistance in the grid cell
        self.neighbors = [] ## the list of neighbors to this grid cell


class agent:
    def __init__(self):
        self.genotype = {"mu_0":0.0,"dev":0.0} ## the dictionary holding important genotypic values
        self.health = 1.0 ## health pool for the agent
        self.resistance = 0.0 ## resistance of the agent to antibacterial
        self.x = 0 ## x coordinate of the agent
        self.y = 0 ## y coordinate of the agent
        self.neighbors = [] ## list of neighbors to this agent
        self.currentCell = None ## the cell this agent currently lives on
    
    def __hash__(self):
        return hash(str(self.x) + "," + str(self.y))
    
    def __eq__(self, other):
        if(hash(self) == hash(other)):
            return True
        return False
    
    ### randomize the genes
    def randomizeGenes(self):
        self.genotype["mu_0"] = np.random.random()
        ## initial draw for resistance
        self.resistance = np.random.normal(loc = self.genotype["mu_0"], scale=self.genotype["dev"])
        if(self.resistance >= 1.0):
            self.resistance = 1.0
        if(self.resistance <= 0.0):
            self.resistance = 0.0
    
    ## place the agent on a cell c
    def place(self, c):
        self.x = c.x
        self.y = c.y
        self.currentCell = c
        c.contents = self
        
    ## The update function updates the agent every iteration
    def update(self, 
               bGrowthRate, 
               nextAgentList, 
               agentsToDie,
               dmgRate,
               dev):
            
        ## Take Damage
        self.damage(dmgRate)
        ## die if health too low
        if(self.health < eps):
            agentsToDie.append(self)
            return
        
        ## degrade concentration
        self.degrade()
        
        ## Try to grow
        offspring = self.grow(bGrowthRate)
        if(offspring != None):
            nextAgentList.append(offspring)
            offspring.mutation()
    
    def mutation(self):
        self.genotype["mu_0"] = self.genotype["mu_0"] + np.random.random()-0.1
        self.genotype["dev"] = self.genotype["dev"] + np.random.random()-0.1
        if(self.genotype["dev"] < 0):
            self.genotype["dev"] = 0.0001
        ## draw for resistance
        self.resistance = np.random.normal(loc = self.genotype["mu_0"], scale=self.genotype["dev"])
        if(self.resistance >= 1.0):
            self.resistance = 1.0
        if(self.resistance <= 0.0):
            self.resistance = 0.0

    
    ## Take damage due to antibacterial in the environment
    def damage(self, dmgRate):
        H_prev = self.health
        c = self.currentCell.concentration
        d = dmgRate
        H_next = H_prev - c*d*H_prev
        self.health = H_next
        
    ## degrade the concentration of antibacterial at this location
    def degrade(self):
        if(self.currentCell.concentration > eps 
           and self.resistance >= eps):
            self.currentCell.concentration -= self.resistance
        elif(self.currentCell.concentration <= eps):
            self.currentCell.concentration = 0.0
    
    ## attempt to grow according to the basal growth rate
    def grow(self, bGrowthRate):
        roll = np.random.random()
        if(roll < bGrowthRate):
            emptyNeighborCells = []
            for c in self.currentCell.neighbors:
                if(c.contents == None):
                    emptyNeighborCells.append(c)
            if(len(emptyNeighborCells) == 0):
                return None
            np.random.shuffle(emptyNeighborCells)
            offspring = agent()
            offspring.place(emptyNeighborCells[0])
            offspring.genotype = self.genotype
            offspring.resistance = np.random.normal(loc = offspring.genotype["mu_0"], scale=offspring.genotype["dev"])
            offspring.genotype["dev"] = self.genotype["dev"]
        else:
            return None
        return offspring

## function to run the sim 
## params: dict of sim parameters
## stats: dict of statistics arrays
def runSim(params, stats):

    ## Set up parameters
    ## Concentration of toxin parameters
    if("initConc" in params):
        initConc = params["initConc"]
    else:
        initConc = 1.0
    if("concPulse" in params):
        concPulse = params["concPulse"]
    else:
        concPulse = 1.0
    if("pulseTime" in params):
        pulseTime = params["pulseTime"]
    else:
        pulseTime = 10


    ## Pop parameters
    if("initPopSize" in params):
        initPopSize = params["initPopSize"]
    else:
        initPopSize = 100

    ## grid parameters
    if("gridWidth" in params):
        gridWidth = params["gridWidth"]
    else:
        gridWidth = 50

    if("gridHeight" in params):
        gridHeight = params["gridHeight"]
    else:
        gridHeight = 50

    ## life parameters
    if("bGrowthRate" in params):
        bGrowthRate = params["bGrowthRate"]
    else:
        bGrowthRate = 0.05
    if("dmgRate" in params):
        dmgRate = params["dmgRate"]
    else:
        dmgRate = 0.1
    if("dev" in params):
        dev = params["dev"]
    else:
        dev = 1.0

    ## sim and stats
    if("numIters" in params):
        numIters = params["numIters"]
    else:
        numIters = 100

    if("pop" in stats):
        pop = stats["pop"]
    else:
        pop = []
    if("births" in stats):
        births = stats["births"]
    else:
        births = []
    if("deaths" in stats):
        deaths = stats["deaths"]
    else:
        deaths = []
    if("avgRes" in stats):
        avgRes = stats["avgRes"]
    else:
        avgRes = []
    if("totalConc" in stats):
        totalConc = stats["totalConc"]
    else:
        totalConc = []
    if("eps" in stats):
        eps = stats["eps"]
    else:
        eps = 10e-8

    ### create grid
    grid = []
    for i in range(gridWidth):
        grid.append([])
        for j in range(gridHeight):
            grid[i].append(cell(i,j))
            grid[i][j].concentration = initConc
    ### Assign neighbors to grid cells
    for i in range(gridWidth):
        for j in range(gridWidth):
            for ii in range(-1,2):
                for jj in range(-1,2):
                    dx = i + ii
                    dy = j + jj
                    if(dx >= 0 and 
                       dx < gridWidth and 
                       dy >= 0 and 
                       dy < gridHeight and
                      ii != 0 and jj != 0):
                        grid[i][j].neighbors.append(grid[dx][dy])
                    
    ### Create initial population
    agentList = []
    for i in range(initPopSize):
        a = agent()
        a.randomizeGenes()
        rx = np.random.randint(0, gridWidth)
        ry = np.random.randint(0, gridHeight)
        while(grid[rx][ry].contents != None):
            rx = np.random.randint(0, gridWidth)
            ry = np.random.randint(0, gridHeight)
        a.place(grid[rx][ry])
        a.randomizeGenes()
        a.genotype["dev"] = dev
        agentList.append(a)

    pulseCounter = pulseTime

    ### Run Sim
    for i in range(numIters):
        
        pulseCounter-=1

        concSum = 0.0
        for x in range(gridWidth):
            for y in range(gridHeight):
                if(pulseCounter <= 0):
                    grid[x][y].concentration += concPulse
                    if(grid[x][y].concentration >= 1.0):
                        grid[x][y].concentration = 1.0
                concSum += grid[x][y].concentration
        totalConc.append(concSum)
        
        resSum = 0.0
        pop.append(len(agentList))
        nextAgentList = [] ### for next iteration population
        agentsToDie = [] ### Agents to remove for next iteration
        np.random.shuffle(agentList)
        for j in range(len(agentList)):
            resSum += agentList[j].resistance
            agentList[j].update(bGrowthRate, 
                                nextAgentList, 
                                agentsToDie,
                                dmgRate,
                                dev)
        births.append(len(nextAgentList))
        deaths.append(len(agentsToDie))

        for j in range(len(agentList)):
            nextAgentList.append(agentList[j])

        agentList = nextAgentList    

        ## remove dead agents
        toDie = set(agentsToDie)
        temp = []
        for a in agentList:
            if a in toDie:
                a.currentCell.contents = None
            else:
                temp.append(a)
        agentList = temp
        
        if(pop[-1] > 0):
            avgRes.append(resSum/float((pop[-1])))
        else:
            avgRes.append(0.0)

        if(pulseCounter <= 0):
            pulseCounter = pulseTime