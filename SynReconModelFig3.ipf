#pragma rtGlobals=3		// Use modern global access method and strict wave access.

CONSTANT tmPnt=19

//Recreates Figure 3 (reduced model)
//Fig 3F is located in wave boundary after running doMany()
//
//To recreate the example graphs for the reduce model, run doOne(6000,1200)
//		The data for Fig 3a is located in the wave titled activity
//		The data for Fig 3b is located in the wave titled predicted
//		The data for Fig 3c is located in the wave titled reservoir
//		The data for Fig 3d is located in the wave titled synapses
function doMany()
	makeConstants()
	make /o/n=(90,90) boundary=NaN
	setscale /p x,0,100,boundary
	setscale /p y,0,100,boundary
	variable synthTm,stimTm
	variable i,j
	for(i=0;i<dimsize(boundary,0);i+=1)
		synthTm=dimoffset(boundary,0)+i*dimdelta(boundary,0)
		for(j=0;j<dimsize(boundary,1);j+=1)
			stimTm=dimoffset(boundary,1)+j*dimdelta(boundary,1)
			doOne(synthTm,stimTm)
			wave predicted
			wavestats /q predicted
			boundary[i][j]=(predicted[numpnts(predicted)-1]-predicted[0])/(v_max-predicted[0])
		endfor
		doUpdate
	endfor
	setscale /p x,dimoffset(boundary,0)/3600,dimdelta(boundary,0)/3600,boundary
	setscale /p y,dimoffset(boundary,1)/3600,dimdelta(boundary,1)/3600,boundary
end

function makeConstants()
	make /o/n=(tmPnt+1) constants
	
	constants[0]=20000		//total number of molecules
	
	constants[1]=10000		//inverse of r1
	constants[2]=1500			//inverse of r3
	
	constants[3]=5e-4			//x0 for r2 sigmoid
	constants[4]=2e-3			//max for r2 sigmoid
	constants[5]=0.057		//x half for r2 sigmoid
	constants[6]=0.003		//rate for r2 sigmoid
	
	constants[7]=5e-4			//x0 for r4 sigmoid
	constants[8]=2e-3			//max for r4 sigmoid
	constants[9]=0.045		//x half for r4 sigmoid
	constants[10]=0.003		//rate for r4 sigmoid
	
	constants[11]=0.01//0.0945		//offset between bound and total synapses
	constants[12]=100			//number of points to calculate slope
	constants[13]=-105		//threshold for reacting to slope
	
//	constants[14]=2.04007959365845	//first increase for activity prediction
//	constants[15]=3.38961720466614	//second increase for activity prediction
	constants[14]=5.4	//first increase for activity prediction
	constants[15]=0	//second increase for activity prediction
	constants[16]=150			//tau for integrating activity
	
	constants[17]=90000//90742		//midpoint for stim and synth
	constants[18]=36000//35642			//beginning point for predictions
	
	constants[tmPnt]=0.1		//time step
end

function doOne(synthTm,stimTm)
	variable synthTm,stimTm
	wave activity,totalN
	wave constants
	wave forModel1
	makeSynth(activity,floor(synthTm/constants[tmPnt]))
	wave newSynth
	variable startPnt=floor(constants[17]-stimTm/constants[tmPnt]/2)
	predictActivity(activity,startPnt,floor(stimTm*0.1))
//	activity=forModel1[p][0]
//	newSynth=forModel1[p][1]
	run(activity,newSynth,totalN)
end

function run(act,syn,totalVal)
	wave act,syn,totalVal
//	variable timer=startMStimer
	wave constants
	duplicate /o totalVal,predicted
	duplicate /o act reservoir,synapses
	STRUCT Molecules molec
	initializeMolecules(molec,totalVal[0])
	reservoir[0]=molec.R
	synapses[0]=molec.S
	wavestats /q totalVal
	variable nMax=v_max
	variable offst=constants[11]
	variable delay=constants[11]/deltax(act)
	variable nNow
	variable pnts=constants[12]
	variable thresh=constants[13]
	variable startNpred=constants[18]
	findLevel /q /EDGE=1 syn,0.5
	variable transition
	if(v_flag)
		transition=numpnts(syn)
	else
		transition=ceil((v_levelX-leftx(syn))/deltax(syn))
	endif
	variable dlt=deltax(act)
	variable endVal
	variable i
	for(i=1;i<numpnts(reservoir);i+=1)
		if(i>startNpred && i<transition)
			nNow=predicted[i-1]
			predicted[i]=predictN(synapses,nMax,nNow,delay,thresh,i-1,dlt,pnts)
		elseif(i==transition)
			endVal=predicted[i-1]
			predicted[i]=endVal
		elseif(i>transition)
			predicted[i]=endVal
		endif
		updateValues(molec,act[i],syn[i],predicted[i])
		synapses[i]=molec.S
		reservoir[i]=molec.R
	endfor
//	print stopMStimer(timer)/1e6
end

function updateValues(molec,act,syn,tots)
	STRUCT Molecules &molec
	variable act,syn,tots
	wave constants
	variable dltS,dltR
	variable r1=1/constants[1]*syn
	variable r3=1/constants[2]
	
	variable r2=constants[4]-constants[3]
	r2/=1+exp((constants[5]-act)/constants[6])
	r2+=constants[3]
	
	variable r4=constants[8]-constants[7]
	r4/=1+exp((constants[9]-act)/constants[10])
	r4+=constants[7]
	
	variable totSyn=tots
	
	variable deltaT=constants[tmPnt]
	
	dltR=-r3*molec.R*(totSyn-molec.S)+r4*molec.S
	dltR*=deltaT
	
	dltS=-(r2+r4)*(molec.S)+r1*molec.F*(totSyn-molec.S)+r3*molec.R*(totSyn-molec.S)
	dltS*=deltaT
	
	molec.R=molec.R+dltR
	molec.S=molec.S+dltS
	molec.F=molec.mx-molec.R-molec.S
end

function predictN(Bwv,nMax,nNow,delay,thresh,val,dlt,pnts)
	wave Bwv
	variable nMax,nNow,delay,thresh,val,dlt,pnts
	variable m=Bwv[val-delay+pnts]-Bwv[val-delay]
	m/=dlt*pnts
	if(m>thresh && m<0)
		m=0
	endif
	return limit(nNow+m*dlt,Bwv[val],nMax)
end

function predictActivity(wv,startPnt,howMany)
	wave wv
	variable ,startPnt,howMany
	wave constants
	wv[constants[18],numpnts(wv)-1]=0
	variable pnts=numpnts(wv)-startPnt+1
	make /o/n=(pnts) predictAct=0
	setscale /p x,0,0.1,predictAct
	variable freq=0.1
	variable up=1/freq/deltax(predictAct)
	variable add
	variable dlt
	variable delta=deltax(predictAct)
	variable val1=constants[14]
	variable val2=constants[15]
	variable startVal=startPnt*deltax(wv)
	variable tau=constants[16]
	variable i
	for(i=1;i<numpnts(predictAct);i+=1)
		if(mod(i-1,up)==0 && floor((i-1)/100)<howMany)
			add=val1
		elseif(mod(i-2,up)==0 && floor((i-2)/100)<howMany)
			add=val2
		else
			add=0
		endif
		dlt=(-predictAct[i-1]+add)*delta/tau
		predictAct[i]=predictAct[i-1]+dlt
	endfor
	setscale /p x,startVal,deltax(wv),predictAct
	wv[startPnt,numpnts(wv)-1]=predictAct[p-startPnt]
end

function makeSynth(wv,howLong)
	wave wv
	variable howLong
	wave constants
	duplicate /o wv newSynth
	newSynth=1
	variable mid=constants[17]
	if(howLong>10)
		newSynth[mid-howLong/2,mid+howLong/2-1]=0
	endif
end

function initializeMolecules(molec,val)
	STRUCT Molecules &molec
	variable val
	wave constants
	molec.S=val
	molec.mx=constants[0]
	molec.R=0
	molec.F=molec.mx-molec.R-molec.S
end

Structure Molecules
	double F
	double S
	double R
	double mx
EndStructure

function getAvgs(wv)
	wave wv
	make /o/n=(dimsize(wv,0),dimsize(wv,2)) one
	setscale /p x,dimoffset(wv,0),dimdelta(wv,0),one
	variable i
	for(i=0;i<dimsize(wv,1);i+=1)
		one=wv[p][i][q]
		getAvg(one)
		wave avg,sd
		make /o/n=(numpnts(avg)/10) holdA,holdS
		setscale /p x,leftx(avg),deltax(avg)*10,holdA,holdS
		holdA=avg(x)
		holdS=sd(x)
		duplicate /o holdA $"avg"+num2str(i)
		duplicate /o holdS $"sd"+num2str(i)
	endfor
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