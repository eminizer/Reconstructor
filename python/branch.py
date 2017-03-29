#Global Variables
MAX_MET_ARRAY_LENGTH = 1
MAX_MU_ARRAY_LENGTH = 50
MAX_EL_ARRAY_LENGTH = 50
MAX_AK4_JET_ARRAY_LENGTH = 100
MAX_AK8_JET_ARRAY_LENGTH = 100
MAX_AK8_SUBJET_ARRAY_LENGTH = 25
MAX_SCALE_ARRAY_LENGTH = 6
MAX_PDF_ARRAY_LENGTH = 100
MAX_ALPHAS_ARRAY_LENGTH = 2

#Imports
from array import array
from ROOT import std

class Branch(object) :

	def __init__(self,readname=None,writename=None,ttreetype='F',inival=-900.,size='1') :
		self.__readname=readname
		self.__writename=writename
		self.__ttreetype=ttreetype
		self.__inival = inival
		if (self.__ttreetype=='I' or self.__ttreetype=='i') and self.__inival==-900. :
			self.__inival=0
		self.__size = size
		#figure out the array type
		self.__arraytype = get_array_type(ttreetype)
		#figure out the array length
		self.__arraylength = get_array_length(size)
		#Set the array to read into
		if readname!=None :
			if ttreetype!='vi' :
				if self.__arraylength>1 :
					self.__readArray = array(self.__arraytype,self.__arraylength*[self.__inival])
				else :
					self.__readArray = array(self.__arraytype,[self.__inival])
			else :
				self.__readArray = std.vector(std.vector('int'))()
		#Set the array to write into
		if writename!=None :
			#If we're just copying over, use the same array to read and write
			if readname!=None :
				self.__writeArray = self.__readArray
			elif self.__arraylength>1 :
				self.__writeArray = array(self.__arraytype,self.__arraylength*[self.__inival])
			else :
				self.__writeArray = array(self.__arraytype,[self.__inival])

	def reset(self) :
		#print 'resetting branch with readname %s and writename %s'%(self.__readname,self.__writename) #DEBUG
		if self.__readname!=None and self.__arraytype!='vi' :
			for i in range(self.__arraylength) :
				self.__readArray[i]=self.__inival
		if self.__writename!=None :
			for i in range(self.__arraylength) :
				self.__writeArray[i]=self.__inival

	def initialize(self,readtree,writetree) :
		self.readtree = readtree
		self.writetree = writetree
		if self.__readname!=None :
			readtree.SetBranchAddress(self.__readname,self.__readArray)
		if self.__writename!=None :
			if self.__arraylength==1 :
				writetree.Branch(self.__writename,self.__writeArray,self.__writename+'/'+self.__ttreetype)
			else :
				writetree.Branch(self.__writename,self.__writeArray,self.__writename+'['+self.__size+']/'+self.__ttreetype)

	def setWriteValue(self,value,index=0) :
		self.__writeArray[index]=value
	def getReadValue(self,index=0) :
		if index<len(self.__readArray) :
			return self.__readArray[index]
		else :
			return self.__inival
	def getReadArray(self) :
		return self.__readArray
	def getWriteValue(self,index=0) :
		return self.__writeArray[index]
	def getReadName(self) :
		return self.__readname
	def getWriteName(self) :
		return self.__writename

def get_array_type(ttreetype) :
	if ttreetype=='vi' :
		return 'vi'
	if ttreetype=='F' :
		return 'f'
	elif ttreetype=='D' :
		return 'd'
	elif ttreetype=='I' :
		return 'i'
	elif ttreetype=='i' :
		return 'I'
	else :
		return None

def get_array_length(size) :
	if size=='1' :
		return 1
	elif size=='met_size' :
		return MAX_MET_ARRAY_LENGTH
	elif size=='mu_size' :
		return MAX_MU_ARRAY_LENGTH
	elif size=='el_size' :
		return MAX_EL_ARRAY_LENGTH
	elif size=='jetAK4CHS_size' :
		return MAX_AK4_JET_ARRAY_LENGTH
	elif size=='jetAK8CHS_size' :
		return MAX_AK8_JET_ARRAY_LENGTH
	elif size=='subjetAK8CHS_size' :
		return MAX_AK8_SUBJET_ARRAY_LENGTH
	elif size=='scale_size' :
		return MAX_SCALE_ARRAY_LENGTH
	elif size=='pdf_size' :
		return MAX_PDF_ARRAY_LENGTH
	elif size=='alphas_size' :
		return MAX_ALPHAS_ARRAY_LENGTH
	else :
		return None