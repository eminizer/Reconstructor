import resource

thismem = {'':0.,'&':0.,'&&':0.,'&&&':0.,'&&&&':0.,'&&&&&':0.}
lastcheckedloc = {'':'INIT','&':'INIT','&&':'INIT','&&&':'INIT','&&&&':'INIT','&&&&&':'INIT'}
mainleaknum = 0
all_changes = {'':[],'&':[],'&&':[],'&&&':[],'&&&&':[],'&&&&&':[]}

class Change(object) :

	def __init__(self,ind,pmem,ploc,cmem,cloc) :
		self._pmem = pmem
		self._ploc = ploc
		self._cmem = cmem
		self._cloc = cloc
		self._memdiff = cmem-pmem
		self._message = ind+' %d change between %s and %s (%d -> %d)'%(self._memdiff,self._ploc,self._cloc,self._pmem,self._cmem)

	def getPMem(self) :
		return self._pmem
	def getCMem(self) :
		return self._cmem
	def printMessage(self) :
		print self._message

def printlevel(ind,startmem=0.,endmem=0.) :
	if not ind in all_changes.keys() :
		return
	for change in all_changes[ind] :
		if change.getPMem()>=startmem and change.getCMem()<=endmem :
			change.printMessage()
			printlevel(ind+'&',change.getPMem(),change.getCMem())

def checkmem(ind,loc,clearlogs=True,) :
	global mainleaknum
	memcheck = float(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
	if memcheck!=thismem[ind] :
		if ind=='' :
			print '----------------'+str(mainleaknum)+'----------------'
			print ind+' %d change between %s and %s (%d -> %d)'%(memcheck-thismem[''],lastcheckedloc[''],loc,thismem[''],memcheck)
			printlevel('&',thismem[ind],memcheck)
			print '----------------'+str(mainleaknum)+'----------------'
			for cl in all_changes.values() :
				while len(cl)>0 :
					cl.pop()
			mainleaknum+=1
		else :
			all_changes[ind].append(Change(ind,thismem[ind],lastcheckedloc[ind],memcheck,loc))
	thismem[ind] = memcheck
	lastcheckedloc[ind] = loc
