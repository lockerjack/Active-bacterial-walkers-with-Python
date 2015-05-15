Peptone=np.zeros([size,size])
for i in range(0,size):
	for j in range(0,size):
		Peptone[i][j] =100*np.exp(-((i-25)**2+(j-25)**2)/169)

for i in range(size):
	for j in range(size):
		Peptone[i,j]=50*(1+sin(2*np.pi*i/(size))*sin(2*np.pi*j/(size)))