import random
f1=open('entrada.in', 'w')
for i in range(10) :
	f1.write('matriz'+str(i)+'\n')
	ord=int(random.uniform(2,+10))
	f2=open('matriz'+str(i), 'w')
	f2.write(str(ord)+'\n')
	for j in range(ord) :
		for k in range(ord) :
			f2.write(str(random.uniform(0,+1000000))+' ')
		f2.write('\n')
	f2.close()
f1.close()
