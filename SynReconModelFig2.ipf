#pragma rtGlobals=3		// Use modern global access method and strict wave access.

CONSTANT zTmPnt=19,dynTmPnt=14
CONSTANT nrnTmPnt=16,numPreCells=2000,connectivity=0.1
CONSTANT samplingTm=0.1

//Recreates Figure 2
//Inputs are the number of repeats (howMany) for each condition
//
//		The data for Fig 2B (purple) is in wave reserv2
//		The data for Fig 2B (yellow) is in wave reserv0
//		The data for Fig 2B (aqua) is in wave reserv1
//		all the error values averages across repeats are located in a wave with the subscript 'sd'
//
//To recreate Fig 2A run runRecon(0,1,0), it is random, but will eventually create a situation whreby the scaffold degades as in the figure
//		The data for Fig 2A top (orange) is located in the wave output[0][ ]
//		The data for Fig 2A bottom (black) is located in the wave output[2][ ]
//		The data for Fig 2A bottom (green) is located in the wave output[7][ ]
function doMany(howMany)
	variable howMany
	variable i
	for(i=0;i<howMany;i+=1)
		runRecon(0,1,0)
		wave avgWeights
		wave output
		duplicate /o avgWeights holdOne
		holdOne=output[8][p]
		duplicate /o holdOne oneReser
		runRecon(0,1,1)
		holdOne=output[8][p]
		concatenate /NP=1 "holdOne;",oneReser
		runRecon(0,0,0)
		holdOne=output[8][p]
		concatenate /NP=1 "holdOne;",oneReser
		if(i==0)
			duplicate /o oneReser manyRecon
		else
			concatenate /NP=2 "oneReser;",manyRecon
		endif
		print time()
	endfor
	reservAvg()
end

function reservAvg()
	wave manyRecon
	make /o/n=(dimsize(manyRecon,0)/100,dimsize(manyRecon,2)) justOne
	setscale /p x,dimoffset(manyRecon,0),dimdelta(manyRecon,0)*100,justOne
	justOne=manyRecon(x)[2][q]
	getAvg(justOne)
	wave avg,sd
	duplicate /o avg reserv0
	duplicate /o sd reserv0sd
	
	justOne=manyRecon(x)[0][q]
	getAvg(justOne)
	duplicate /o avg reserv1
	duplicate /o sd reserv1sd
	
	justOne=manyRecon(x)[1][q]
	getAvg(justOne)
	duplicate /o avg reserv2
	duplicate /o sd reserv2sd
end

function getAvg(wv)
	wave wv
	make /o/n=(dimsize(wv,0)) avg,sd
	setscale /p x,dimoffset(wv,0),dimdelta(wv,0),avg,sd
	make /o/n=(dimsize(wv,1)) hold
	variable i
	for(i=0;i<numpnts(avg);i+=1)
		hold=wv[i][p]
		wavestats /q hold
		avg[i]=v_avg
		sd[i]=v_sdev/sqrt(v_npnts)
	endfor
end

function runRecon(syn,midStim,long)
	variable syn,midStim,long
//	variable timer=startMStimer
	
	variable numCells=200//binomialNoise(numPreCells,connectivity)
	
	makeNeuronConstants()
	makeZconstants()
	makeSTDPminVisConstants()
	makeDynamicsConstants(numCells)
	wave nrnConst,zConst,STDPconst,dynConst
	variable delta=nrnConst[nrnTmPnt]
	struct NeuronParams neuron0
	initializeNeuron(neuron0)
	
	make /o/n=(3,2,1) postSdet
	make /o/n=(3,1,numCells) preSdet
	initializeDetectors(preSdet,postSdet)
	make /o/n=(numCells) deltaWeights,weights=zConst[16],allZ
	
	make /o/n=(numCells,7,2) synapses
	initializeSynapses(synapses,weights)
	NVAR dopa,synthesis

	initializeProteins(synapses)
	wave proteinStatus,turnOver,whatCond
	
	make /o/n=(numCells) spikeTms,whichCell
	make /o/n=0 outputSpikes
	
	make /o/n=(0) avgWeights,totalN,totalB,totalK4,totalLoss,totalGain,activity
	make /o/n=(dimsize(synapses,1)+4,0) output
	make /o/n=(numCells) synHolder
	make /o/n=(dimsize(synapses,1)) oneSynTime
	make /o/n=(5,1) activities
	make /o/n=(dimsize(activities,0)) oneActSet
	
	variable freq=0.1
	variable stim20=1200
	variable stimLong=3600
	variable tmAfter=1000
	
	variable pnts20=stim20*freq
	variable pntsLong=stimLong*freq
	variable pntsAfter=tmAfter*freq

	
	variable lastSpikePnt=0
	variable startTm=0
	lastSpikePnt=runLowFreq(spikeTms,whichCell,outputSpikes,synapses,weights,deltaWeights,preSdet,postSdet,neuron0,delta,freq,lastSpikePnt,startTm,pnts20,numCells,1,proteinStatus,turnOver,whatCond)

	startTm=ceil(neuron0.T+zConst[18]*4)
	lastSpikePnt=runHighFreq(spikeTms,whichCell,outputSpikes,synapses,weights,deltaWeights,preSdet,postSdet,neuron0,delta,100,100,startTm,lastSpikePnt,proteinStatus,turnOver,whatCond)
	
	dopa=1
	
	startTm=ceil(neuron0.T+zConst[18]*4)
	lastSpikePnt=runLowFreq(spikeTms,whichCell,outputSpikes,synapses,weights,deltaWeights,preSdet,postSdet,neuron0,delta,freq,lastSpikePnt,startTm,ceil(60*freq),numCells,1,proteinStatus,turnOver,whatCond)
	
	dopa=0
	
	startTm=ceil(neuron0.T+zConst[18]*4)
	lastSpikePnt=runLowFreq(spikeTms,whichCell,outputSpikes,synapses,weights,deltaWeights,preSdet,postSdet,neuron0,delta,freq,lastSpikePnt,startTm,pnts20,numCells,1,proteinStatus,turnOver,whatCond)
	
	if(long)
		startTm=ceil(neuron0.T+zConst[18]*4)
		lastSpikePnt=runLowFreq(spikeTms,whichCell,outputSpikes,synapses,weights,deltaWeights,preSdet,postSdet,neuron0,delta,freq,lastSpikePnt,startTm,pntsLong,numCells,1,proteinStatus,turnOver,whatCond)
	else
		advanceLIF(outputSpikes,synapses,weights,deltaWeights,preSdet,postSdet,neuron0,3600/delta,delta,proteinStatus,turnOver,whatCond)
	endif
	
	if(syn==0)
		synthesis=0
	endif
	
	if(midStim)
		if(long==0)
			lastSpikePnt=neuron0.T/delta
		endif
	
		startTm=8464
		lastSpikePnt=runLowFreq(spikeTms,whichCell,outputSpikes,synapses,weights,deltaWeights,preSdet,postSdet,neuron0,delta,freq,lastSpikePnt,startTm,pnts20,numCells,1,proteinStatus,turnOver,whatCond)
	else
		advanceLIF(outputSpikes,synapses,weights,deltaWeights,preSdet,postSdet,neuron0,3600/delta,delta,proteinStatus,turnOver,whatCond)
	endif
	
	advanceLIF(outputSpikes,synapses,weights,deltaWeights,preSdet,postSdet,neuron0,2400/delta,delta,proteinStatus,turnOver,whatCond)
	
	if(syn==0)
		synthesis=1
	endif
	
//	advanceLIF(outputSpikes,synapses,weights,deltaWeights,preSdet,postSdet,neuron0,100/delta,delta,proteinStatus,turnOver,whatCond)
	
	lastSpikePnt=neuron0.T/delta
	startTm=18000
	lastSpikePnt=runLowFreq(spikeTms,whichCell,outputSpikes,synapses,weights,deltaWeights,preSdet,postSdet,neuron0,delta,freq,lastSpikePnt,startTm,pntsAfter,numCells,1,proteinStatus,turnOver,whatCond)
	
//	advanceLIF(outputSpikes,synapses,weights,deltaWeights,preSdet,postSdet,neuron0,10/delta,delta,proteinStatus,turnOver,whatCond)
	
	setscale /p x,0,samplingTm/3600,avgWeights,totalN,totalB,totalK4,totalLoss,totalGain,activity
	setscale /p y,0,samplingTm/3600,output,activities
	make /o/n=2 line1,line2
	setscale /p x,dimoffset(output,1),dimoffset(output,1)+(dimsize(output,1)-1)*dimdelta(output,1),line1,line2
	line1=1
	line2=-1
		
//	print stopMStimer(timer)/1e6/neuron0.T
//	print neuron0.T
end

function makeZconstants()
	make /o/n=(zTmPnt+1) zConst
	zConst[0]=1.3*0.25		//a x meta
	zConst[1]=0.95*0.25		//a y meta
	zConst[2]=3.5*0.25		//a y tilt
	zConst[3]=3.5*0.25		//a z tilt
	
	zConst[4]=sqrt(1e-4)		//sigma x
	zConst[5]=sqrt(1e-4)		//sigma y
	zConst[6]=sqrt(1e-4)		//sigma z
	
	zConst[7]=exp(-1)			//theta for gyz
	
	zConst[8]=200				//tau x
	zConst[9]=200				//tau y
	zConst[10]=200			//tau z
	zConst[11]=600			//tau gamma for gyx
	
	zConst[12]=1				//tau up for PRP
	zConst[13]=1000			//tau down for PRP
	
	zConst[14]=1				//constant for STDP weight change
	
	zConst[15]=0.33			//fraction of synapses intialized to up state
	
	zConst[16]=0.05			//value of weight if x is -1
	zConst[17]=3				//ratio of max weight to w0

	zConst[18]=0.003			//sd of spike timing with stimulation
	
	zConst[zTmPnt]=0.1		//time step
	
	zConst[4]*=zConst[8]/sqrt(zConst[zTmPnt])
	zConst[5]*=zConst[9]/sqrt(zConst[zTmPnt])
	zConst[6]*=zConst[10]/sqrt(zConst[zTmPnt])
end

function makeDynamicsConstants(numCells)
	variable numCells
	make /o/n=(dynTmPnt+1) dynConst
	dynConst[0]=0						//Threshold to define when z is up

	dynConst[1]=10000				//time to go from free to bound (s/number free)
	dynConst[2]=1500					//time to go from reserve to bound (s/number reserve)
	
	dynConst[3]=5e-4					//x0 of sigmoid for activity function 1
	dynConst[4]=2e-3					//max of sigmoid for activity function 1
	dynConst[5]=0.057					//x half of sigmoid for activity function 1
	dynConst[6]=0.003					//rate of sigmoid for activity function 1
	
	dynConst[7]=5e-4					//x0 of sigmoid for activity function 2
	dynConst[8]=2e-3					//max of sigmoid for activity function 2
	dynConst[9]=0.045					//x half of sigmoid for activity function 2
	dynConst[10]=0.003				//rate of sigmoid for activity function 2
	
	dynConst[11]=150					//decay time for activity
	
	dynConst[12]=numCells*100		//total number of molecules
	
	dynConst[13]=0					//starting reservoir value
	
	dynConst[dynTmPnt]=0.1			//time step
end

function makeNeuronConstants()
	make /o/n=(nrnTmPnt+1) nrnConst
	nrnConst[0]=-0.07				//Reversal potential (V)
	nrnConst[1]=-0.05				//Threshold potential (V)
	nrnConst[2]=-0.07				//Reset potential (V)
	nrnConst[3]=0.1				//Threshold increase after a spike (V)
	nrnConst[4]=0.002				//Refractory threshold decay (s)
	nrnConst[5]=0.02				//Membrane time constant for excitatory cell (s)
	nrnConst[6]=0.01				//Membrane time constant for inhibitory cell (s)
	nrnConst[7]=0					//Excitatory conductance reversal potential (V)
	nrnConst[8]=0.005				//AMPA conductance decay time (s)
	nrnConst[9]=0.1				//NMDA conductance decay time (s)
	nrnConst[10]=0.5				//AMPA contribution to excitatory conductance
	nrnConst[11]=1-nrnConst[10]	//NMDA contribution to excitatory conductance
	nrnConst[12]=-0.08			//Inhibitory conductance reversal potential (V)
	nrnConst[13]=0.01				//GABA conductance decay time (s)
	nrnConst[14]=10				//Post spike adaptive magnitude (?)
	nrnConst[15]=0.25				//Post spike adaptive decay (s)
	
	nrnConst[nrnTmPnt]=1e-4		//Time step
end

function saveInfo(weightsWv,synWv,protWv,tmsWv)
	wave weightsWv,synWv,protWv,tmsWv
	wave avgWeights,output,totalN,allZ,totalB,totalK4,totalGain,totalLoss,whatCond,activity
//	wave activities,oneActSet
	NVAR reservoir,synthesis,pool,gainR,lossR
	make /o/n=(numpnts(avgWeights)+1) avgWeights
	make /o/n=(numpnts(totalN)+1) totalN
	make /o/n=(numpnts(totalB)+1) totalB
	make /o/n=(numpnts(totalB)+1) totalK4
	make /o/n=(numpnts(totalGain)+1) totalGain
	make /o/n=(numpnts(totalLoss)+1) totalLoss
	make /o/n=(numpnts(activity)+1) activity
	wavestats /q/m=1 weightsWv
	avgWeights[numpnts(avgWeights)-1]=v_avg
	allZ=whatCond[p]||protWv[p]
	totalN[numpnts(totalN)-1]=sum(allZ)
	totalB[numpnts(totalB)-1]=sum(protWv)
	allZ=tmsWv[p][1]*protWv[p]
	totalK4[numpnts(totalK4)-1]=sum(allZ)
	totalGain[numpnts(totalGain)-1]=gainR
	totalLoss[numpnts(totalLoss)-1]=lossR
	allZ=synWv[p][6][0]
	wavestats /q/m=1 allZ
	activity[numpnts(activity)-1]=v_avg
	make /o/n=(dimsize(output,0)) oneSynTime
	if(dimsize(output,1)==0)
		oneSynTime=NaN
//		oneActSet=NaN
	else
		oneSynTime[0,dimsize(synWv,1)-1]=synWv[0][p][0]
		oneSynTime[numpnts(oneSynTime)-4]=protWv[0][0]
		oneSynTime[numpnts(oneSynTime)-3]=reservoir
		oneSynTime[numpnts(oneSynTime)-2]=synthesis
		oneSynTime[numpnts(oneSynTime)-1]=pool
//		oneActSet=synWv[p][6][0]
	endif
	concatenate /NP=1 "oneSynTime;",output
//	concatenate /NP=1 "oneActSet;",activities
end

function runLowFreq(spkTms,whchCll,outputSpks,synWv,weightsWv,dltWeightWv,preWv,postWv,nrn,dlt,freq,lastSpikePnt,strtTm,pnts,numCells,reps,proteinWv,tmsWv,whichWv)
	wave spkTms,whchCll,outputSpks,synWv,weightsWv,dltWeightWv,preWv,postWv,proteinWv,tmsWv,whichWv
	Struct NeuronParams &nrn
	variable dlt,freq,lastSpikePnt,strtTm,pnts,numCells,reps
	variable nextPnt
	variable i,j,k
	for(i=1;i<=pnts;i+=1)
		for(k=0;k<reps;k+=1)
		makeLowFreqSpikes(spkTms,whchCll,i/freq+strtTm+k*0.05)
			for(j=0;j<numCells;j+=1)
				nextPnt=round(spkTms[j]/dlt)-lastSpikePnt
				lastSpikePnt=round(spkTms[j]/dlt)
				advanceLIF(outputSpks,synWv,weightsWv,dltWeightWv,preWv,postWv,nrn,nextPnt,dlt,proteinWv,tmsWv,whichWv)
				advanceSynapseAfterSpike(synWv,preWv,postWv,weightsWv,dltWeightWv,tmsWv,nrn,whchCll[j])
			endfor
		endfor
	endfor
	return lastSpikePnt
end

function runHighFreq(spkTms,whchCll,outputSpks,synWv,weightsWv,dltWeightWv,preWv,postWv,nrn,dlt,freq,howMany,strtTm,lastSpikePnt,proteinWv,tmsWv,whichWv)
	wave spkTms,whchCll,,outputSpks,synWv,weightsWv,dltWeightWv,preWv,postWv,proteinWv,tmsWv,whichWv
	Struct NeuronParams &nrn
	variable dlt,freq,howMany,strtTm,lastSpikePnt
	variable nextPnt
	makeHighFreqSpikes(spkTms,whchCll,strtTm,freq,howMany)
	wave allSpikes,allCells
	variable i
	for(i=0;i<numpnts(allSpikes);i+=1)
		nextPnt=round(allSpikes[i]/dlt)-lastSpikePnt
		lastSpikePnt=round(allSpikes[i]/dlt)
		advanceLIF(outputSpks,synWv,weightsWv,dltWeightWv,preWv,postWv,nrn,nextPnt,dlt,proteinWv,tmsWv,whichWv)
		advanceSynapseAfterSpike(synWv,preWV,postWv,weightsWv,dltWeightWv,tmsWv,nrn,allCells[i])
	endfor
	return lastSpikePnt
end

function advanceLIF(outputWv,synWv,weightsWv,dltWeightWv,preWv,postWv,nrn,howMany,dlt,protWv,tmsWv,whatWv)
	wave outputWv,synWv,weightsWv,dltWeightWv,preWv,postWv,protWv,tmsWv,whatWv
	struct NeuronParams &nrn
	variable howMany,dlt
	wave zConst
	variable i
	for(i=0;i<howMany;i+=1)
		updatePotential(nrn,0)
		if(mod(nrn.T,zConst[zTmPnt])<dlt)
			updateProteinPresence(synWv,protWv,tmsWv,whatWv)
			updateSynapse(synWv,protWv,tmsWv)
			updateWeights(weightsWv,synWv)
			updateDecaysTms(synWv,tmsWv,-1)
			if(mod(nrn.T,samplingTm)<zConst[zTmPnt])
				saveInfo(weightsWv,synWv,protWv,tmsWv)
			endif
		endif
		if(nrn.S)
			nrn.S=0
			make /o/n=(numpnts(outputWv)+1) $NameOfWave(outputWv)
			outputWv[numpnts(outputWv)-1]=nrn.T
			getDeltaWeight(preWv,postWv,dltWeightWv,1,nrn.T)
			updateDetectors(postWv,-1,nrn.T)
			updateSynapseAfterPostSpike(synWv,dltWeightWv,tmsWv)
			updateWeights(weightsWv,synWv)
		endif
	endfor
end

function advanceSynapseAfterSpike(synWv,preWv,postWv,weightsWv,dltWeightWv,tmsWv,nrn,which)
	wave synWv,preWv,postWv,weightsWv,dltWeightWv,tmsWv
	struct NeuronParams &nrn
	variable which
	preExcSpikeHappens(nrn,weightsWv[which])
	variable dltW=getDeltaWeight(preWv,postWv,dltWeightWv,0,nrn.T)
	updateDetectors(preWv,which,nrn.T)
	updateSynapseAfterPreSpike(synWv,tmsWv,which,dltW)
	updateWeights(weightswv,synWv)
end

function updateWeights(weightwv,synWv)
	wave weightWv,synWv
	wave zConst
	weightWv=synWv[p][0]*(zConst[17]-1)
	weightWv+=1+zConst[17]
	weightWv*=zConst[16]*0.5
end

function updateDecaysTms(synWv,timesWv,val)
	wave synWv,timesWv
	variable val
	wave dynConst
	if(val>0)
		timesWv[val][0]=dynConst[4]-dynConst[3]
		timesWv[val][0]/=1+exp((dynConst[5]-synWv[p][6][0])/dynConst[6])
		timesWv[val][0]+=dynConst[3]
		
		timesWv[val][1]=dynConst[8]-dynConst[7]
		timesWv[val][1]/=1+exp((dynConst[9]-synWv[p][6][0])/dynConst[10])
		timesWv[val][1]+=dynConst[7]
	else
		timesWv[][0]=dynConst[4]-dynConst[3]
		timesWv[][0]/=1+exp((dynConst[5]-synWv[p][6][0])/dynConst[6])
		timesWv[][0]+=dynConst[3]
		
		timesWv[][1]=dynConst[8]-dynConst[7]
		timesWv[][1]/=1+exp((dynConst[9]-synWv[p][6][0])/dynConst[10])
		timesWv[][1]+=dynConst[7]
	endif
end

function updateSynapse(synWv,protWv,timesWv)
	wave synWv,protWv,timesWv
	wave zConst
	wave dynConst
	updateXvals(synWv)
	updateYvals(synWv)
	updateZvals(synWv,protWv)
	updateGyx(synWv)
	updateGzy(synWv)
	updateActivity(synWv)
	updateDecaysTms(synWv,timesWv,-1)
	synWv[][0,3][1]/=zConst[q+8]
	synWv[][][1]*=zConst[zTmPnt]
	synWv[][][0]+=synWv[p][q][1]
	synWv[][4][0]=(synWv[p][3][0]>=zConst[7])
end

function updateSynapseAfterPreSpike(synWv,timesWv,which,dltW)
	wave synWv,timesWv
	variable which,dltW
	wave zConst
	wave dynConst
	wave nrnConst
	synWv[which][0][1]=limit((synWv[p][0][0]-synWv[p][2][0]),0,inf)
	synWv[which][0][1]*=zConst[14]
	synWv[which][0][1]+=1
	synWv[which][0][1]*=-dltW
	synWv[which][0][1]*=(1+synWv[p][0][0])
	
	synWv[which][3][1]=(synWv[p][2][0]>=synWv[p][0][0])
	synWv[which][3][1]*=dltW
	synWv[which][3][1]*=(1-synWv[p][3][0])
		
	synWv[which][0][0]+=synWv[p][0][1]
	synWv[which][3][0]+=synWv[p][3][1]
	synWv[which][4][0]=(synWv[p][3][0]>=zConst[7])
	
	synWv[which][6][0]+=dltW
	updateDecaysTms(synWv,timesWv,which)
end

function updateSynapseAfterPostSpike(synWv,deltaWv,timesWv)
	wave synWv,deltaWv,timesWv
	wave zConst
	wave dynConst
	wave nrnConst
	synWv[][0][1]=limit((synWv[p][2][0]-synWv[p][0][0]),0,inf)
	synWv[][0][1]*=zConst[14]
	synWv[][0][1]+=1
	synWv[][0][1]*=deltaWv[p]
	synWv[][0][1]*=(1-synWv[p][0][0])
	
	synWv[][3][1]=(synWv[p][0][0]>=synWv[p][2][0])
	synWv[][3][1]*=deltaWv[p]
	synWv[][3][1]*=(1-synWv[p][3][0])
		
	synWv[][0][0]+=synWv[p][0][1]
	synWv[][3][0]+=synWv[p][3][1]
	synWv[][4][0]=(synWv[p][3][0]>=zConst[7])
	
	synWv[][6][0]+=deltaWv[p]
	updateDecaysTms(synWv,timesWv,-1)
end

function updateXvals(synWv)
	wave synWv
	wave zConst
	synWv[][0][1]=-synWv[p][0][0]*(synWv[p][0][0]-1)*(synWv[p][0][0]+1)
	synWv[][0][1]+=gnoise(1)*zConst[4]
	synWv[][0][1]+=zConst[0]*(1-synWv[p][4][0])*(synWv[p][1][0]-synWv[p][0][0])
end

function updateYvals(synWv)
	wave synWv
	wave zConst
	synWv[][1][1]=-synWv[p][1][0]*(synWv[p][1][0]-1)*(synWv[p][1][0]+1)
	synWv[][1][1]+=gnoise(1)*zConst[5]
	synWv[][1][1]+=zConst[2]*(synWv[p][4][0])*(synWv[p][0][0]-synWv[p][1][0])
	synWv[][1][1]+=zConst[1]*(1-synWv[p][5][0])*(synWv[p][2][0]-synWv[p][1][0])
end

function updateZvals(synWv,protWv)
	wave synWv,protWv
	wave zConst
	synWv[][2][1]=-synWv[p][2][0]*(synWv[p][2][0]-1)*(synWv[p][2][0]+1)
	synWv[][2][1]+=gnoise(1)*zConst[6]
	synWv[][2][1]+=zConst[3]*(synWv[p][5][0])*(synWv[p][1][0]-synWv[p][2][0])
	
	synWv[][2][1]+=(protWv[p]-1)*(synWv[p][2][0]>0)
end

function updateGyx(synWv)
	wave synWv
	synWv[][3][1]=-synWv[p][3][0]
end

function updateGzy(synWv)
	wave synWv
	wave zConst
	NVAR dopa
	synWv[][5][1]=dopa/zConst[12]*(1-synWv[p][5][0])
	synWv[][5][1]-=(synWv[p][5][0])/zConst[13]
end

function updateActivity(synWv)
	wave synWv
	wave dynConst
	synWv[][6][1]=-synWv[p][6][0]/dynConst[11]
end

function initializeSynapses(synWv,weightsWv)
	wave synWv,weightsWv
	wave zConst
	wave dynConst
	synWv[][3,6][0]=0
	make /o/n=(dimsize(synWv,0)) order=p,random
	random=gnoise(1)
	sort random,random,order
	variable val=numpnts(order)*zConst[15]
	variable i
	for(i=0;i<numpnts(order);i+=1)
		if(i<val-1)
			synWv[order[i]][0,2][0]=1
		else
			synWv[order[i]][0,2][0]=-1
		endif
	endfor
	
	synWv[0][0,2][0]=-1
	
	updateWeights(weightswv,synWv)
	variable /g dopa=0
	variable /g synthesis=1
	variable /g reservoir=dynConst[13]
	variable /g gainR=0,lossR=0
	variable /g pool=dynConst[12]-dynConst[13]
	killwaves /z order,random
end

function makeLowFreqSpikes(spikeTmsWv,whichCellWv,tm)
	wave spikeTmsWv,whichCellWv
	variable tm
	wave zConst
	whichCellWv=p
	spikeTmsWv=gnoise(zConst[18])+tm
	sort spikeTmsWv,spikeTmsWv,whichCellWv
end

function makeHighFreqSpikes(spikeTmsWv,whichCellWv,startTm,freq,howMany)
	wave spikeTmsWv,whichCellWv
	variable startTm,freq,howMany
	wave zConst
	variable i
	for(i=0;i<howMany;i+=1)
		whichCellWv=p
		spikeTmsWv=gnoise(zConst[18])+startTm+i*(1/freq)
		if(i==0)
			duplicate /o spikeTmsWv allSpikes
			duplicate /o whichCellWv allCells
		else
			concatenate /NP=0 NameOfWave(spikeTmsWv)+";",allSpikes
			concatenate /NP=0 NameOfWave(whichCellWv)+";",allCells
		endif
	endfor
	sort allSpikes,allSpikes,allCells
	killwaves /z other
end

function updateProteinPresence(synWv,protWv,turnOverWv,whatWv)
	wave synWv,protWv,turnOverWv,whatWv
	wave dynConst
	variable delta=dynConst[dynTmPnt]
	make /o/w/u/n=(numpnts(protWv)) order=p
	whatWv=synWv[p][2][0]>dynConst[0]
	variable val
	NVAR reservoir
	NVAR synthesis
	NVAR pool
	NVAR gainR
	gainR=0
	NVAR lossR
	lossR=0
	variable probF,probR
	variable i,j
	for(i=0;i<numpnts(protWv);i+=1)
		j=getRandI(numpnts(order))
		val=order[j]
		if(protWv[val]==0 && whatWv[val]==1)
			probF=(synthesis*pool/dynConst[1]*delta)>abs(enoise(1))
			probR=(reservoir/dynConst[2]*delta)>abs(enoise(1))
			if(probF&&probR)
				protWv[val]=1
				if(enoise(1)>=0)
					reservoir-=1
					lossR+=1
				else
					pool-=1
				endif
			else
				protWv[val]=probF+probR
				reservoir-=probR
				lossR+=probR
				pool-=probF
			endif
		elseif(protWv[val])
			probF=(turnOverWv[val][0]*delta)>abs(enoise(1))
			probR=(turnOverWv[val][1]*delta)>abs(enoise(1))
			if(probF&&probR)
				protWv[val]=0
				if(enoise(1)>=0)
					reservoir+=1
					gainR+=1
				else
					pool+=1
				endif
			else
				protWv[val]=!(probF+probR)
				reservoir+=probR
				gainR+=probR
				pool+=probF
			endif
		endif
		deletepoints j,1,order
	endfor
end

function initializeProteins(synWv)
	wave synWv
	NVAR pool
	wave dynConst
	variable numCells=dimsize(synWv,0)
	make /o/b/u/n=(numCells) proteinStatus=0, whatCond=0
	proteinStatus[]=(synWv[p][2][0]==1)
	
	pool-=sum(proteinStatus)
			
	make /o/n=(numCells,2) turnOver
	turnOver[][0]=dynConst[3]
	turnOver[][1]=dynConst[7]
end

function getRandI(val)
	variable val
	return limit(floor(abs(enoise(val))),0,val-1)
end

function placeNLs(wv,which)
	wave wv
	variable which
	wave dynConst
	make /o/n=(dimsize(wv,1)) hold
	setscale /p x,dimoffset(wv,1),dimdelta(wv,1),hold
	hold=wv[which][p]
	wavestats /q hold
	variable val
	val=x2pnt(hold,v_maxLoc-.1)
	make /o/n=(val) hold1
	hold1=wv[which][p]
	make /o/n=100 hist1,hist2
	wavestats /q hold1
	setscale /i x,0,v_avg+3*v_sdev,hist1
	histogram /b=2 hold1 hist1
	wavestats /q hist1
	hist1/=v_sum
	wavestats /q hold
	val=x2pnt(hold,v_maxLoc+.1)
	make /o/n=(numpnts(hold)-val) hold1
	hold1=hold[p+val]
	wavestats /q hold1
	setscale /i x,0,v_avg+3*v_sdev,hist2
	histogram /b=2 hold1 hist2
	wavestats /q hist2
	hist2/=v_sum
	
	duplicate /o hist2 NL1,NL2
	NL1=dynConst[4]-dynConst[3]
	NL1/=1+exp((dynConst[5]-x)/dynConst[6])
	NL1+=dynConst[3]
	
	NL2=dynConst[8]-dynConst[7]
	NL2/=1+exp((dynConst[9]-x)/dynConst[10])
	NL2+=dynConst[7]
	killwaves /z hold1
end

function updatePotential(nrn,extInp)
	struct NeuronParams &nrn
	variable extInp
	wave nrnConst
	updateConductances(nrn)
	variable leak=(nrnConst[0]-nrn.V)
	variable synaptic=(nrnConst[7]-nrn.V)*(nrn.Aval*nrnConst[10]+nrn.Nval*nrnConst[11])
	synaptic+=(nrnConst[12]-nrn.V)*(nrn.Adval+nrn.Gval)
	variable delta=(leak+synaptic+extInp)*nrnConst[nrnTmPnt]/nrnConst[5]
	nrn.V+=delta
	variable threshold=nrnConst[1]+nrn.Rval
	if(nrn.V>=threshold)
		nrn.V=nrnConst[2]
		nrn.S=1
		postSpikeHappens(nrn)
	endif
	nrn.T+=nrnConst[nrnTmPnt]
end

function updateConductances(nrn)
	struct NeuronParams &nrn
	nrn.Aval*=nrn.Aupdate
	nrn.Nval+=(nrn.Aval-nrn.Nval)*nrn.Nupdate
	nrn.Gval*=nrn.Gupdate
	nrn.Adval*=nrn.Adupdate
	nrn.Rval*=nrn.Rupdate
end

function preExcSpikeHappens(nrn,wgt)
	struct NeuronParams &nrn
	variable wgt
	nrn.Aval+=wgt
end

function preInhSpikeHappens(nrn,wgt)
	struct NeuronParams &nrn
	variable wgt
	nrn.Gval+=wgt
end

function postSpikeHappens(nrn)
	struct NeuronParams &nrn
	wave nrnConst
	nrn.Rval=nrnConst[3]
	nrn.Adval+=nrnConst[14]
end

function initializeNeuron(nrn)
	struct NeuronParams &nrn
	wave nrnConst
	nrn.V=nrnConst[0]
	nrn.S=0
	nrn.Rval=0
	nrn.Rupdate=exp(-nrnConst[nrnTmPnt]/nrnConst[4])
	nrn.Aval=0
	nrn.Aupdate=exp(-nrnConst[nrnTmPnt]/nrnConst[8])
	nrn.Nval=0
	nrn.Nupdate=nrnConst[nrnTmPnt]/nrnConst[9]
	nrn.Gval=0
	nrn.Gupdate=exp(-nrnConst[nrnTmPnt]/nrnConst[13])
	nrn.Adval=0
	nrn.Adupdate=exp(-nrnConst[nrnTmPnt]/nrnConst[15])
	
	nrn.T=-nrnConst[nrnTmPnt]
end

structure NeuronParams
	float V				//current membrane potential (mV)
	char S				//Spike (1=Yes, 0=No)
	float Rval			//Refractory value
	variable Rupdate	//Refractory update value
	float Aval			//AMPA conductance value
	variable Aupdate	//AMPA update multiplier
	float Nval			//NMDA conductance value
	variable Nupdate	//NMDA update multiplier
	float Gval			//GABA conductance value
	variable Gupdate	//GABA update multiplier
	float Adval			//Adaptation conductance value
	variable Adupdate	//Adaptation update multiplier
	
	variable T			//time counter
endStructure

function makeSTDPminHipConstants()
	make /o/n=(8) STDPconst
	STDPconst[0]=0.0168	//tau plus
	STDPconst[1]=0		//tau x
	STDPconst[2]=0.0337	//tau minus
	STDPconst[3]=0.04		//tau y
	STDPconst[4]=3.5e-3	//A minus 2
	STDPconst[5]=0		//A minus 3
	STDPconst[6]=5.3e-3	//A plus 2
	STDPconst[7]=8e-3	//A plus 3
end

function makeSTDPminVisConstants()
	make /o/n=(8) STDPconst
	STDPconst[0]=0.0168	//tau plus
	STDPconst[1]=0		//tau x
	STDPconst[2]=0.0337	//tau minus
	STDPconst[3]=0.114//0.04//	//tau y
	STDPconst[4]=2e-4//7.1e-3	//A minus 2
	STDPconst[5]=0		//A minus 3
	STDPconst[6]=0		//A plus 2
	STDPconst[7]=5e-4//6.5e-3	//A plus 3
end

function makeSTDPminLorConstants()
	make /o/n=(8) STDPconst
	STDPconst[0]=0.0168	//tau plus
	STDPconst[1]=0		//tau x
	STDPconst[2]=0.0337	//tau minus
	STDPconst[3]=0.04//	//tau y
	STDPconst[4]=2e-4//7.1e-3	//A minus 2
	STDPconst[5]=0		//A minus 3
	STDPconst[6]=0		//A plus 2
	STDPconst[7]=5e-4//6.5e-3	//A plus 3
end

function getDeltaWeight(preWv,postWv,deltaWeightsWv,preOrpost,tm)
	wave preWv,postWv,deltaWeightsWv
	variable preOrpost,tm
	wave STDPconst
	variable val
	variable dltW
	if(preOrpost==0)		//presynaptic spike
		val=getDetectorVal(preWv,postWv,preOrpost,tm)
		dltW=val*(STDPconst[4])
		return dltW
	elseif(preOrPost==1)			//postsynaptic spike
		val=getDetectorVal(preWv,postWv,preOrPost,tm)
		wave preDetVals
		deltaWeightsWv=preDetVals[p]*(STDPconst[6]+STDPconst[7]*val)
	endif
end

function getDetectorVal(preWv,postWv,whichDet,tm)
	wave preWv,postWv
	variable whichDet,tm
	variable tau,x0,A
	variable val
	if(whichDet)			//need to get all presynaptic values
		wave preDetVals
		preDetVals=(preWv[2][0][p]*exp(-(tm-preWv[1][0][p])/preWv[0][0][p]))/preWv[0][0][p]
	endif
	tau=postWv[0][whichDet]
	x0=postWv[1][whichDet]
	A=postWv[2][whichDet]
	val=(A*exp(-(tm-x0)/tau))/tau
	return val
end

function updateDetectors(detValsWv,whichCell,tm)
	wave detValsWv
	variable whichCell,tm
	if(whichCell<0)			//postsynaptic spike
		updateIndivDetectors(detValsWv,0,0,tm)
		updateIndivDetectors(detValsWv,1,0,tm)
	else							//presynaptic spike
		updateIndivDetectors(detValsWv,0,whichCell,tm)
	endif
end

function updateIndivDetectors(detValsWv,whichDet,whichCell,when)
	wave detValsWv
	variable whichDet,whichCell,when
	variable tau=detValsWv[0][whichDet][whichCell]
	variable x0=detValsWv[1][whichDet][whichCell]
	variable A=detValsWv[2][whichDet][whichCell]
	detValsWv[1][whichDet][whichCell]=when
	detValsWv[2][whichDet][whichCell]=1+A*exp(-(when-x0)/tau)
end

function initializeDetectors(preWv,postWv)
	wave preWv,postWv
	wave STDPconst
	preWv=0
	postWv=0
	preWv[0][0][]=STDPconst[0]
	postWv[0][0][]=STDPconst[2]
	postWv[0][1][]=STDPconst[3]
	
	make /o/n=(dimsize(preWv,2)) preDetVals
end