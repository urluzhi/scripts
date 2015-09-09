#RNA input 
import sys
RNA_seq=sys.argv[1]
RNA_seq=RNA_seq.upper().replace('T','U')

keylist=['CG','GC','AU','UA','GU','UG']
valuelist=[3,3,2,2,2,2]
D=dict(zip(keylist,valuelist))

trace=list('.'*len(RNA_seq))

#return the maximum pairing in the case of two branches
def biggest(i,j):
	multi={k:hair(i,k)+hair(k+1,j) for k in range(i+1,j)}
	for k in multi.keys():
		if multi[k] == max(multi.values()):
			return [k,multi[k]]


#j>i in default
#only consider the simplest case: hairpin
def hair(i,j):
	pair=RNA_seq[i-1]+RNA_seq[j-1]
	if j-i<=3:
		return 0
	else:
		if pair in D.keys():
			# pairs[i]=j
			return D[pair]+hair(i+1,j-1)
		if not pair in D.keys():
			return fill(i+1,j-1)

def fill(i,j):
	pair=RNA_seq[i-1]+RNA_seq[j-1]
	if j-i<=3:
		return 0
	else:
		#if just add max([fill(i,k)+fill(k+1,j) for k in range(i+1,j)]) will cost too much time 
		#since every pair will consider the possibility of forming two branches
		if pair in D.keys():
			if D[pair]+fill(i+1,j-1) > biggest(i,j)[1]:
				return D[pair]+fill(i+1,j-1)
			else:
				return biggest(i,j)[1]
				
		if not pair in D.keys():
			return max(fill(i+1,j),fill(i,j-1),biggest(i,j)[1])

def track(i,j):
	pair=RNA_seq[i-1]+RNA_seq[j-1]
	if j-i>3 and pair in D.keys():
		if fill(i,j)==D[pair]+fill(i+1,j-1):
			trace[i-1]='('
			trace[j-1]=')'
			track(i+1,j-1)
		if fill(i,j)==biggest(i,j)[1]:
			track(i,biggest(i,j)[0])
			track(biggest(i,j)[0]+1,j)			
	elif j-i>3 and not pair in D.keys():
			if fill(i,j)==fill(i+1,j):
				track(i+1,j)				
			if fill(i,j)==fill(i,j-1):
				track(i,j-1)				
			if fill(i,j)==biggest(i,j)[1]:
				track(i,biggest(i,j)[0])
				track(biggest(i,j)[0]+1,j)
				
print sys.argv[0]
print "sequence: "+RNA_seq
print "length: "+str(len(RNA_seq))
print RNA_seq
track(1,len(RNA_seq))
print ''.join(trace)
print "max H-bonds: "+str(fill(1,len(RNA_seq)))	
