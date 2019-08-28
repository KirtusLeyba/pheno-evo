package main

import pes "./PhenoEvoSim"
import "flag"
import "time"

func main() {

	tempSeed := time.Now().UnixNano() / int64(time.Millisecond)

	//Set Global Variables from parameters
	epsilonParam := flag.Float64("epsilon", 10e-8, "For small value float comparisons.")
	basalGrowthRateParam := flag.Float64("basalGrowthRate", 1.0, "The base probability for agent growth.")
	envNoiseParam := flag.Float64("envNoise", 0.2, "The maximum noise with which to perturb the toxin signal.")
	tradeOffParam := flag.Float64("tradeoff", 1.0, "The tradeoff coefficient for resistance cost.")
	alphaParam := flag.Float64("alpha", 1.0, "Parameterizes the cost function of increased resistance.")
	toxinConcParam := flag.Float64("toxinConc", 1.0, "The amount of toxin administered per pulse.")
	diluteRateParam := flag.Int("diluteRate", 100, "currentPop/(diluteRate) = (number of agents to keep per dilution)")
	widthParam := flag.Int("width", 50, "The width in patches of the simulation space.")
	heightParam := flag.Int("height", 50, "The height in patches of the simulation space.")
	initPopSizeParam := flag.Int("initPopSize", 200, "The initial population size of agents.")
	switchRateParam := flag.Float64("switchRate", 0.4, "The phenotypic switching probability of the initial population.")
	global_w0Param := flag.Float64("w0", 0.1, "The weight value of the 1st Gamma Distribution.")
	global_m0Param := flag.Float64("m0", 0.2, "The mean parameter of the 1st Gamma Distribution.")
	global_m1Param := flag.Float64("m1", 0.8, "The mean parameter of the 2nd Gamma Distribution.")
	global_v0Param := flag.Float64("v0", 0.01, "The variance parameter of the 1st Gamma Distribution.")
	global_v1Param := flag.Float64("v1", 0.01, "The variance parameter of the 2nd Gamma Distribution.")
	envResponseParam := flag.Float64("envResponse", 0.5, "The probability for an agent to switch its degrade rate to match the signal at its location.")
	pulseTicksParam := flag.Int("pulseTicks", 20, "The number of ticks between toxin pulses.")
	seedParam := flag.Int64("seed", tempSeed, "The seed for the random number generation.")
	itersParam := flag.Int("iters", 1000, "The number of iterations to run the simulation.")
	pulseMinXParam := flag.Int("pulseMinX", 21, "The minimum X point to pulse from.")
	pulseMaxXParam := flag.Int("pulseMaxX", 29, "The maximum X point to pulse from.")
	pulseMinYParam := flag.Int("pulseMinY", 0, "The minimum Y point to pulse from.")
	pulseMaxYParam := flag.Int("pulseMaxY", 49, "The maximum Y point to pulse from.")

	flag.Parse()

	pes.Epsilon = *(epsilonParam)
	pes.BasalGrowthRate = *(basalGrowthRateParam)
	pes.EnvNoise = *(envResponseParam)
	pes.TradeOff = *(tradeOffParam)
	pes.Alpha = *(alphaParam)
	pes.ToxinConc = *(toxinConcParam)
	pes.DiluteRate = *(diluteRateParam)
	pes.Width = *(widthParam)
	pes.Height = *(heightParam)
	pes.InitPopSize = *(initPopSizeParam)
	pes.SwitchRate = *(switchRateParam)
	pes.Global_W0 = *(global_w0Param)
	pes.Global_M0 = *(global_m0Param)
	pes.Global_V0 = *(global_v0Param)
	pes.Global_M1 = *(global_m1Param)
	pes.Global_V1 = *(global_v1Param)
	pes.EnvResponse = *(envResponseParam)
	pes.PulseTicks = *(pulseTicksParam)
	pes.Seed  = *(seedParam)
	pes.EnvNoise = *(envNoiseParam)
	pes.PulseMinX = *(pulseMinXParam)
	pes.PulseMinY = *(pulseMinYParam)
	pes.PulseMaxX = *(pulseMaxXParam)
	pes.PulseMaxY = *(pulseMaxYParam)

	if(pes.PulseMaxY >= pes.Height) {
		pes.PulseMaxY = pes.Height-1
	}
	if(pes.PulseMaxX >= pes.Width) {
		pes.PulseMaxY = pes.Width-1
	}

	numIters := *(itersParam)

	pes.RunSim(numIters)
}