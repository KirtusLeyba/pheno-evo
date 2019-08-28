package PhenoEvoSim

import "gonum.org/v1/gonum/stat/distuv"
import "golang.org/x/exp/rand"
import mr "math/rand"
import "math"
import "fmt"


//Global Variables
var Epsilon float64 //for small number comparisons
var BasalGrowthRate float64
var EnvNoise float64
var TradeOff float64
var Alpha float64
var ToxinConc float64
var DiluteRate int
var GenoVar float64
var Width int
var Height int
var InitPopSize int
var SwitchRate float64
var Global_W0 float64
var Global_M0 float64
var Global_V0 float64
var Global_M1 float64
var Global_V1 float64
var EnvResponse float64
var PulseTicks int
var Seed int64
var PulseMinX int
var PulseMinY int
var PulseMaxX int
var PulseMaxY int

type Agent struct {
	X, Y int
	M_0, V_0, W_0, M_1, V_1 float64 //distribution variables
	DegradeRate float64
	SwitchRate float64
	EnvResponse float64
	Health float64
	PatchHere *Patch
	GrowthRate float64
}

type Patch struct {
	X, Y int
	Toxin float64
	Signal float64
	AgentHere *Agent
}


// draw a random variable from the defined gamma distribution
func GammaDraw(alpha float64, beta float64, seed uint64) float64 {
	source := rand.NewSource(seed)
	g := distuv.Gamma{alpha, beta, source}
	return g.Rand()
}

// Draw a random phenotype from a mixed gamma distribution
func DrawPhenotype(w_0 float64, m_0 float64, v_0 float64, m_1 float64, v_1 float64) float64 {
	result := 0.0
	seed := mr.Uint64()
	if mr.Float64() < w_0 {
		result = GammaDraw(m_0*m_0/v_0, 1/(v_0/m_0), seed) 
	} else {
		result = GammaDraw(m_1*m_1/v_1, 1/(v_1/m_1), seed)
	}
	if(result >= 1.0) {
		result = 1.0
	}
	if(result <= 0.0) {
		result = 0.0
	}
	return result
}

// apply poison to the grid within a defined region
func PoisonGrid(xMin, xMax, yMin, yMax int, toxinConc float64, grid [][]*Patch) {
	for i:=xMin; i <= xMax; i++ {
		for j:=yMin; j <= yMax; j++ {
			patch := grid[i][j]
			patch.Toxin += toxinConc
		}
	}
}

// diffuse toxin through the grid
func DiffuseGrid(grid [][]*Patch, width int, height int, diffRate float64) {

	//need to go through patches in a random order
	xIndeces := []int{}
	yIndeces := []int{}
	for i :=0; i < width; i++ {
		xIndeces = append(xIndeces, i)
	}
	for i := 0; i < height; i++ {
		yIndeces = append(yIndeces, i)
	}
	//shuffle indeces
	mr.Shuffle(len(xIndeces), func(i,j int){xIndeces[i], xIndeces[j] = xIndeces[j], xIndeces[i]})
	mr.Shuffle(len(yIndeces), func(i,j int){yIndeces[i], yIndeces[j] = yIndeces[j], yIndeces[i]})

	for ii:=0; ii < width; ii++ {
		for jj:=0; jj < height; jj++ {
			i := xIndeces[ii]
			j := yIndeces[jj]
			sourcePatch := grid[i][j]
			sourceToxin := sourcePatch.Toxin
			targetPatches := []*Patch{}


			for dx:=-1;dx<=1;dx++ {
				for dy:=-1;dy<=1;dy++{
					x := i + dx
					y := j + dy
					if(x >= 0 && x < width && y >= 0 && y < height) {
						if(!(x == 0 && y == 0)) {
							targetPatches = append(targetPatches, grid[x][y])
						}						
					}
				}
			}
			mr.Shuffle(len(targetPatches), func(iz, jz int) { targetPatches[iz], targetPatches[jz] = targetPatches[jz], targetPatches[iz] })
			for k := 0; k < len(targetPatches); k++ {
				dT := (sourceToxin/9.0)*diffRate
				targetPatch := targetPatches[k]
				targetPatch.Toxin += dT
				sourcePatch.Toxin -= dT
			}
		}
	}
}

func SignalPatches(grid [][]*Patch, width int, height int) {
	for i:=0; i<width; i++ {
		for j:=0; j<height; j++ {
			patch := grid[i][j]
			patch.Signal = patch.Toxin + (mr.Float64()*EnvNoise) - 0.5*EnvNoise
		}
	}
}

//dilute the simulation environment
// changes agentList in place
func DiluteGrid(agentList *([]*Agent), grid *([][]*Patch)) {
	newPopNum := int(math.Floor(float64(len(*agentList))/float64(DiluteRate)))
	mr.Shuffle(len(*agentList), func(i, j int) { (*agentList)[i], (*agentList)[j] = (*agentList)[j], (*agentList)[i] })
	newAgentList := (*agentList)[0:newPopNum]
	for i:=newPopNum; i < len(*(agentList)); i++ {
		(*grid)[(*agentList)[i].X][(*agentList)[i].Y].AgentHere = nil
	}
	*agentList = newAgentList
}


func InitializeSim() ([][](*Patch), [](*Agent)) {

	mr.Seed(Seed)

	// create the simulation grid
	grid := [][](*Patch){}
	for i:=0;i<Width;i++ {
		grid = append(grid, [](*Patch){})
		for j:=0;j<Height;j++ {
			newPatch := Patch{}
			newPatch.X = i
			newPatch.Y = j
			newPatch.Toxin = 0.0
			newPatch.Signal = 0.0
			grid[i] = append(grid[i], &newPatch)
		}
	}

	// create the initial population of agents
	agentList := [](*Agent){}
	for i:=0;i<InitPopSize;i++{
		newAgent := Agent{}
		rx := mr.Intn(Width)
		ry := mr.Intn(Height)
		for(grid[rx][ry].AgentHere != nil) {
			rx = mr.Intn(Width)
			ry = mr.Intn(Height)
		}
		newAgent.X = rx
		newAgent.Y = ry
		newAgent.Health = 1.0
		newAgent.SwitchRate = SwitchRate
		newAgent.EnvResponse = EnvResponse
		newAgent.W_0 = Global_W0
		newAgent.M_0 = Global_M0
		newAgent.V_0 = Global_V0
		newAgent.M_1 = Global_M1
		newAgent.V_1 = Global_V1
		newAgent.PatchHere = grid[rx][ry]
		newAgent.DegradeRate = DrawPhenotype(newAgent.W_0, newAgent.M_0, newAgent.V_0, newAgent.M_1, newAgent.V_1)
		newAgent.GrowthRate = BasalGrowthRate * (1.0 - TradeOff * math.Pow(newAgent.DegradeRate, Alpha))
		agentList = append(agentList, &newAgent)

	}

	return grid, agentList

}

func RunSim(iterations int) {
	// first setup the environment
	grid, agentList := InitializeSim()

	header := "iter,x,y,agentHere,Toxin,agentDegradeRate\n"
	fmt.Printf(header)

	for iter:=0; iter < iterations; iter++ {

		//signal the environment
		SignalPatches(grid, Width, Height)

		//pulse the environment
		if(iter%PulseTicks == 0) {
			PoisonGrid(PulseMinX,PulseMaxX,PulseMinY,PulseMaxY, ToxinConc, grid)
		}


		for i:=0;i<Width;i++ {
			for j:=0;j<Height;j++{
				agentHere := 0
				toxinHere := grid[i][j].Toxin
				agentDegradeRateHere := 0.0
				if(grid[i][j].AgentHere != nil){
					agentHere = 1
					agentDegradeRateHere = grid[i][j].AgentHere.DegradeRate
				}
				//print output
				fmt.Printf("%d,%d,%d,%d,%f,%f\n",iter, i, j, agentHere, toxinHere, agentDegradeRateHere)
			}
		}

		//Switch phenotypes for agents
		for i:=0; i < len(agentList); i++ {
			agent := agentList[i]
			if mr.Float64() < agent.SwitchRate {
				if(mr.Float64() < agent.EnvResponse) {
					agent.DegradeRate = agent.PatchHere.Signal
				} else {
					agent.DegradeRate = DrawPhenotype(agent.W_0, agent.M_0, agent.V_0, agent.M_1, agent.V_1)
				}
			}
		}

		// effect of toxin on agents and vice versa
		temp := []*Agent{}
		for i:=0; i < len(agentList); i++ {

			agent := agentList[i]
			patchHere := agent.PatchHere
			if(patchHere.Toxin >= Epsilon) {
				agent.Health -= patchHere.Toxin
			} else {
				agent.Health = 1.0 //recover health immediately if no toxin
			}
			if(agent.Health >= Epsilon) {
				temp = append(temp, agent)
				if(patchHere.Toxin > Epsilon) {
					patchHere.Toxin -= agent.DegradeRate
				}
				if(patchHere.Toxin <= Epsilon) {
					patchHere.Toxin = 0.0
				}
			} else {
				patchHere.AgentHere = nil
			}
		}

		//reproduction
		for i:=0; i < len(agentList); i++ {
			growthSpace := []*Patch{}
			agent := agentList[i]
			if(mr.Float64() < agent.GrowthRate){
				x := agent.X
				y := agent.Y

				for dx:=-1;dx<=1;dx++ {
					for dy:=-1;dy<=1;dy++{
						xx := x + dx
						yy := y + dy
						if(!(dx == 0 && dy == 0)){
							if(xx >= 0 && xx < Width && yy >= 0 && yy < Height) {
								if(grid[xx][yy].AgentHere == nil){
									growthSpace = append(growthSpace, grid[xx][yy])
								}
							}
						}
					}
				}
				if(len(growthSpace) > 0) {
					mr.Shuffle(len(growthSpace), func(i, j int) { growthSpace[i], growthSpace[j] = growthSpace[j], growthSpace[i] })
					nx := growthSpace[0].X
					ny := growthSpace[0].Y

					//create offspring
					Offspring := Agent{}
					Offspring.X = nx
					Offspring.Y = ny
					Offspring.Health = 1.0
					Offspring.SwitchRate = agent.SwitchRate
					Offspring.EnvResponse = agent.EnvResponse
					Offspring.DegradeRate = agent.DegradeRate
					Offspring.W_0 = agent.W_0
					Offspring.M_0 = agent.M_0
					Offspring.V_0 = agent.V_0
					Offspring.M_1 = agent.M_1
					Offspring.V_1 = agent.V_1
					Offspring.PatchHere = grid[nx][ny]
					Offspring.DegradeRate = DrawPhenotype(Offspring.W_0, Offspring.M_0, Offspring.V_0, Offspring.M_1, Offspring.V_1)
					Offspring.GrowthRate = BasalGrowthRate * (1.0 - TradeOff * math.Pow(Offspring.DegradeRate, Alpha))
					grid[nx][ny].AgentHere = &Offspring
					temp = append(temp, &Offspring)

				}
			}




		}

		//diffuse toxin
		// shuffle list to remove update bias
		DiffuseGrid(grid, Width, Height, 0.5)

		//replace old population
		agentList = temp

		//dilute if neccessary
		if(len(agentList) >= Width*Height) {
			DiluteGrid(&agentList, &grid)
		}

	}

}