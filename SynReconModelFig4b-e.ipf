#pragma rtGlobals=3		// Use modern global access method and strict wave access.

//Recreates Figures 4c - f
//Inputs are the number of repeats (howMany) for each figure
//
//		The data for Fig 4b are located in waves recon3 (with PSI) and recon0 (without PSI)
//		The data for Fig 4c are located in waves recon1 (with PSI) and recon4 (without PSI)
//		The data for Fig 4d are located in waves recon2 (with PSI) and recon5 (without PSI)
//		The data for Fig 4e are located in waves full1 (with PSI) and full0 (without PSI)
//		the error values averaged across repeats are located in a wave titled 'reconsd#,' where the # correspond to the number at the end of the wave of interest
function doMany(howMany)
	variable howMany
	variable i
	for(i=0;i<howMany;i+=1)
		doRecon()
		wave oneRecon,oneFull
		if(i==0)
			duplicate /o oneRecon allRecon
			duplicate /o oneFull allFull
		else
			concatenate /NP=3 "oneRecon;",allRecon
			concatenate /NP=3 "oneFull;",allFull
		endif
		print time()
	endfor
	doAvgs()
end

function doRecon()
	makeConstants()
	wave constants
	make /o/n=(600*constants[tmPnt]) aTime1,bTime1,pTime1,dTime1,aTime2,bTime2,pTime2,dTime2,cTime,act1r2,act1r4,act2r2,act2r4
	setscale /p x,0,1/constants[tmPnt],aTime1,bTime1,pTime1,dTime1,aTime2,bTime2,pTime2,dTime2,cTime,act1r4,act2r2,act2r4
	
	wave act,noAct,actLong,actCont
	act2r2=0
	act2r4=0
	
	wave synthRecon
	duplicate /o synthRecon synth
	
	runOne(noAct)
	wave avgWeights1,avgWeights2,reserv
	concatenate /NP=1 "avgWeights2;reserv;",avgWeights1
	duplicate /o avgWeights1 oneRecon
	
	runOne(act)
	concatenate /NP=1 "avgWeights2;reserv;",avgWeights1
	concatenate /NP=2 "avgWeights1;",oneRecon
	
	runOne(actLong)
	concatenate /NP=1 "avgWeights2;reserv;",avgWeights1
	concatenate /NP=2 "avgWeights1;",oneRecon
	
	synth=1
	
	runOne(noAct)
	concatenate /NP=1 "avgWeights2;reserv;",avgWeights1
	concatenate /NP=2 "avgWeights1;",oneRecon
	
	runOne(act)
	concatenate /NP=1 "avgWeights2;reserv;",avgWeights1
	concatenate /NP=2 "avgWeights1;",oneRecon
	
	runOne(actLong)
	concatenate /NP=1 "avgWeights2;reserv;",avgWeights1
	concatenate /NP=2 "avgWeights1;",oneRecon
	
	make /o/n=(600*constants[tmPnt]) aTime1,bTime1,pTime1,dTime1,aTime2,bTime2,pTime2,dTime2,cTime,act1r2,act1r4,act2r2,act2r4
	setscale /p x,0,1/constants[tmPnt],aTime1,bTime1,pTime1,dTime1,aTime2,bTime2,pTime2,dTime2,cTime,act1r4,act2r2,act2r4

	wave synthContFull
	duplicate /o synthContFull synth
	
	runOne(actCont)
	concatenate /NP=1 "avgWeights2;reserv;",avgWeights1
	duplicate /o avgWeights1 oneFull
	
	synth=1
	runOne(actCont)
	concatenate /NP=1 "avgWeights2;reserv;",avgWeights1
	concatenate /NP=2 "avgWeights1;",oneFull
end

function doAvgs()
	wave allRecon,allFull
	make /o/n=(dimsize(allRecon,0),dimsize(allRecon,3)) all
	setscale /p x,0,1/3600,all
	variable i
	for(i=0;i<dimsize(allRecon,2);i+=1)
		all=allRecon[p][0][i][q]
		getAvg(all)
		wave avg,sd
		make /o/n=(numpnts(avg)/10) hold
		setscale /p x,leftx(avg),deltax(avg)*10,hold
		hold=avg(x)
		duplicate /o hold $"recon"+num2str(i)
		hold=sd(x)
		duplicate /o hold $"reconSd"+num2str(i)
	endfor
	
	make /o/n=(dimsize(allFull,0),dimsize(allFull,3)) all
	setscale /p x,0,1/3600,all
	for(i=0;i<dimsize(allFull,2);i+=1)
		all=allFull[p][0][i][q]
		getAvg(all)
		wave avg,sd
		make /o/n=(numpnts(avg)/10) hold
		setscale /p x,leftx(avg),deltax(avg)*10,hold
		hold=avg(x)
		duplicate /o hold $"full"+num2str(i)
		hold=sd(x)
		duplicate /o hold $"fullSd"+num2str(i)
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

function runOne(thisAct)
	wave thisAct
	wave aTime1,aTime2,bTime1,bTime2,cTime,dTime1,dTime2,pTime1,pTime2,act1R2,act1r4
	make /o/n=0 no
	makeAtime(aTime1,{20})
	makeAtime(aTime2,no)
	makeBtime(bTime1,no)
	makeBtime(bTime2,no)
	makeCtime(cTime,{20})
	makePDtime(dTime1,no)
	makePDtime(dTime2,no)
	makePDtime(pTime1,{20})
	makePDtime(pTime2,no)
	makeActs(act1r2,act1r4,thisAct)
	runModel()
end

function makeConstants()
	make /o/n=(tmPnt+1) constants
	
	constants[0]=0.8				//starting fraction of synapses in weak basal state
	
	constants[1]=1/60				//basal rate for changing weak to strong (inv min) (alpha)
	constants[2]=1/15				//basal rate for changing strong to weak (inv min) (beta)
	constants[3]=10				//beta level following LFS
	constants[4]=4					//duration of elevation of beta following LFS (min)
	constants[5]=50				//tm1 for duration of d or p following LFS or HFS
	constants[6]=10				//tm2 for duration of d or p following LFS or HFS
	constants[7]=30				//tm1 for duration of c following sLFS or sHFS
	constants[8]=30				//tm2 for duration of c following sLFS or sHFS
	
	constants[9]=1/60				//tau early
	constants[10]=1e-4			//tau late bound
	constants[11]=1/10			//tau late unbound
	
	constants[12]=1/3000			//r1
	constants[13]=1/30			//r3
	
	constants[14]=100				//multiplyer between # synapses and total pool size
	
	constants[15]=0.02			//x0 of sigmoid for activity function 1
	constants[16]=0.25				//max of sigmoid for activity function 1
	constants[17]=0.058			//x half of sigmoid for activity function 1
	constants[18]=0.003			//rate of sigmoid for activity function 1
	
	constants[19]=0.02			//x0 of sigmoid for activity function 2
	constants[20]=0.25				//max of sigmoid for activity function 2
	constants[21]=0.045			//x half of sigmoid for activity function 2
	constants[22]=0.003			//rate of sigmoid for activity function 2
	
	constants[tmPnt]=60			//1 divided by the sampling time in min
end

function runModel()
	wave constants
	wave aTime1,bTime1,pTime1,dTime1,act1r2,act1r4
	wave aTime2,bTime2,pTime2,dTime2,act2r2,act2r4
	wave cTime,synth
	variable tauE=constants[9]
	variable tauL=constants[10]
	variable tauLf=constants[11]
	variable r1=constants[12]
	variable r3=constants[13]
	variable r1val,r3val
	
	duplicate /o aTime1 avgWeights1,reserv
	duplicate /o aTime2 avgWeights2
	avgWeights1=0
	avgWeights2=0
	
	make /o/n=(1000) synapses1,synapses2,weak,strong
	variable val1=round(numpnts(synapses1)*constants[0])
	synapses1[0,val1-1]=2
	synapses1[val1,numpnts(synapses1)-1]=3
	variable start1=val1+(numpnts(synapses1)-val1)*2
	
	variable val2=round(numpnts(synapses2)*constants[0])
	synapses2[0,val2-1]=2
	synapses2[val2,numpnts(synapses2)-1]=3
	variable start2=val2+(numpnts(synapses2)-val2)*2
	
	variable /g pool=(numpnts(synapses1)+numpnts(synapses2))*constants[14]
	variable /g reservoir=0
	
	variable delta=deltax(avgWeights1)
	variable prob1,prob2
	
	variable i
	for(i=0;i<numpnts(avgWeights1);i+=1)
		r1val=r1*pool*delta*synth[i]
		r3val=r3*reservoir*delta
		updateSynapses(synapses1,aTime1[i]*delta,bTime1[i]*delta,cTime[i]*delta,dTime1[i]*delta,pTime1[i]*delta,tauE*delta,tauL*delta,tauLf*delta,r1val,act1r2[i]*delta,r3val,act1r4[i]*delta)
		weak=synapses1[p]<3
		strong=synapses1[p]>=3
		avgWeights1[i]=1/start1*(sum(weak)+2*sum(strong))*100
		
		updateSynapses(synapses2,aTime2[i]*delta,bTime2[i]*delta,cTime[i]*delta,dTime2[i]*delta,pTime2[i]*delta,tauE*delta,tauL*delta,tauLf*delta,r1val,act2r2[i]*delta,r3val,act2r4[i]*delta)
		weak=synapses2[p]<3
		strong=synapses2[p]>=3
		avgWeights2[i]=1/start2*(sum(weak)+2*sum(strong))*100
		
		reserv[i]=reservoir
	endfor
end

function makeAtime(wv,tms)
	wave wv,tms
	wave constants
	wv=constants[1]
	variable i
	for(i=0;i<numpnts(tms);i+=1)
		wv[x2pnt(wv,tms[i])]+=constants[tmPnt]
	endfor
end

function makeBtime(wv,tms)
	wave wv,tms
	wave constants
	wv=constants[2]
	variable i
	for(i=0;i<numpnts(tms);i+=1)
		wv[x2pnt(wv,tms[i]),x2pnt(wv,tms[i]+constants[4])]=constants[3]
	endfor
end

function makePDtime(wv,tms)
	wave wv,tms
	wave constants
	wv=0
	variable i
	for(i=0;i<numpnts(tms);i+=1)
		wv+=limit((x-tms[i])/constants[5]*exp(1-(x-tms[i])/constants[6]),0,inf)
	endfor
end

function makeCtime(wv,tms)
	wave wv,tms
	wave constants
	wv=0
	variable i
	for(i=0;i<numpnts(tms);i+=1)
		wv+=limit((x-tms[i])/constants[7]*exp(1-(x-tms[i])/constants[8]),0,inf)
	endfor
end

function makeActs(wv1,wv2,actWv)
	wave wv1,wv2,actWv
	wave constants
	wv1=constants[16]-constants[15]
	wv1/=1+exp((constants[17]-actWv[p])/constants[18])
	wv1+=constants[15]
	
	wv2=constants[20]-constants[19]
	wv2/=1+exp((constants[21]-actWv[p])/constants[22])
	wv2+=constants[19]
end

function updateSynapses(wv,aVal,bVal,cVal,dVal,pVal,tauE,tauL,tauLf,r1,r2,r3,r4)
	wave wv
	variable aVal,bVal,cVal,dVal,pVal,tauE,tauL,tauLf,r1,r2,r3,r4
	NVAR pool,reservoir
	variable prob1,prob2,prob3
	variable thisOne
	duplicate /o wv order
	order=p
	variable val
	variable i,j
	for(i=0;i<numpnts(wv);i+=1)
		j=getRandI(numpnts(order))
		val=order[j]
		switch(wv[val])
			case 0:
				if(tauL>abs(enoise(1)))
					wv[val]=2
				endif	
				break
			case 1:
				prob1=tauE>abs(enoise(1))
				prob2=cVal>abs(enoise(1))
				if(prob1 && prob2)
					if(enoise(1)>0)
						wv[val]=2
					else
						wv[val]=0
					endif
				elseif(prob1)
					wv[val]=2
				elseif(prob2)
					wv[val]=0
				endif
				break
			case 2:
				prob1=dVal>abs(enoise(1))
				prob2=aVal>abs(enoise(1))
				if(prob1 && prob2)
					if(enoise(1)>0)
						wv[val]=3
					else
						wv[val]=1
					endif
				elseif(prob1)
					wv[val]=1
				elseif(prob2)
					wv[val]=3
				endif
				break
			case 3:
				prob1=bVal>abs(enoise(1))
				prob2=pVal>abs(enoise(1))
				if(prob1 && prob2)
					if(enoise(1)>0)
						wv[val]=2
					else
						wv[val]=4
					endif
				elseif(prob1)
					wv[val]=2
				elseif(prob2)
					wv[val]=4
				endif
				break
			case 4:
				prob1=tauE>abs(enoise(1))
				prob2=cVal>abs(enoise(1))
				if(prob1 && prob2)
					if(enoise(1)>0)
						wv[val]=3
					else
						wv[val]=5
					endif
				elseif(prob1)
					wv[val]=3
				elseif(prob2)
					wv[val]=5
				endif
				break
			case 5:
				prob1=tauLf>abs(enoise(1))
				prob2=r1>abs(enoise(1))
				prob3=r3>abs(enoise(1))
				if(prob1 && prob2 && prob3)
					thisOne=abs(enoise(3))
					if(thisOne<=1)
						wv[val]=3
					elseif(thisOne<=2)
						wv[val]=6
						reservoir-=1
					else
						wv[val]=6
						pool-=1
					endif
				elseif(prob1 && prob2)
					if(enoise(1)>0)
						wv[val]=3
					else
						pool-=1
						wv[val]=6
					endif
				elseif((prob1 && prob3))
					if(enoise(1)>0)
						reservoir-=1
						wv[val]=6
					else
						wv[val]=3
					endif
				elseif(prob2 && prob3)
					if(enoise(1)>0)
						pool-=1
						wv[val]=6
					else
						reservoir-=1
						wv[val]=6
					endif
				elseif(prob1)
					wv[val]=3
				elseif(prob2)
					pool-=1
					wv[val]=6
				elseif(prob3)
					reservoir-=1
					wv[val]=6
				endif
				break
			case 6:
				prob1=tauL>abs(enoise(1))
				prob2=r2>abs(enoise(1))
				prob3=r4>abs(enoise(1))
				if(prob1 && prob2 && prob3)
					thisOne=abs(enoise(3))
					if(thisOne<=1)
						wv[val]=3
						if(enoise(1)>0)
							reservoir+=1
						else
							pool+=1
						endif
					elseif(thisOne<=2)
						wv[val]=5
						reservoir+=1
					else
						wv[val]=5
						pool+=1
					endif
				elseif(prob1 && prob2)
					if(enoise(1)>0)
						wv[val]=3
						if(enoise(1)>0)
							reservoir+=1
						else
							pool+=1
						endif
					else
						pool+=1
						wv[val]=5
					endif
				elseif((prob1 && prob3))
					if(enoise(1)>0)
						reservoir+=1
						wv[val]=5
					else
						wv[val]=3
						if(enoise(1)>0)
							reservoir+=1
						else
							pool+=1
						endif
					endif
				elseif(prob2 && prob3)
					if(enoise(1)>0)
						pool+=1
						wv[val]=5
					else
						reservoir+=1
						wv[val]=5
					endif
				elseif(prob1)
					wv[val]=3
					if(enoise(1)>0)
						reservoir+=1
					else
						pool+=1
					endif
				elseif(prob2)
					pool+=1
					wv[val]=5
				elseif(prob3)
					reservoir+=1
					wv[val]=5
				endif
				break
		endSwitch
		deletepoints j,1,order
	endfor
end

function getRandI(val)
	variable val
	return limit(floor(abs(enoise(val))),0,val-1)
end