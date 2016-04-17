#pragma rtGlobals=3		// Use modern global access method and strict wave access.

//To recreate Fig5a switch the below code with the corresponding code from the Fig 1 file
function updateXvals(synWv)
	wave synWv
	wave zConst
	synWv[][0][1]=-synWv[p][0][0]*(synWv[p][0][0]-1)*(synWv[p][0][0]+1)
	synWv[][0][1]+=gnoise(1)*zConst[4]*(synWv[p][6][0]>0.01)
	synWv[][0][1]+=zConst[0]*(1-synWv[p][4][0])*(synWv[p][1][0]-synWv[p][0][0])
end

function updateYvals(synWv)
	wave synWv
	wave zConst
	synWv[][1][1]=-synWv[p][1][0]*(synWv[p][1][0]-1)*(synWv[p][1][0]+1)
	synWv[][1][1]+=gnoise(1)*zConst[5]*(synWv[p][6][0]>0.01)
	synWv[][1][1]+=zConst[2]*(synWv[p][4][0])*(synWv[p][0][0]-synWv[p][1][0])
	synWv[][1][1]+=zConst[1]*(1-synWv[p][5][0])*(synWv[p][2][0]-synWv[p][1][0])
end

//To recreate Fig5b switch the below code with the corresponding code from the Fig 1 file
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
				prob2=(aVal>abs(enoise(1)))*(r4>(0.025*0.0166667))
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
				prob1=(bVal>abs(enoise(1)))*(r4>(0.025*0.0166667))
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