import threading
import time


def calc_squares(n):
	print('Calculando os quadrados')
	for i in n:
		time.sleep(0.2)
		print('Quadrado: ',i*i)
		
		
def calc_cubes(n):
	print('Calculando os cubos')
	for i in n:
		time.sleep(0.2)
		print('Cubo: ',i*i*i)		
		

n = [1,2,3,4]


t = time.time()

# ~ calc_squares(n)
# ~ calc_cubes(n)

t1 = threading.Thread(target = calc_squares,args = (n,))
t2 = threading.Thread(target = calc_cubes,args = (n,))

t1.start()
t2.start()

# ~ t1.join()
# ~ t2.join()
		
print(time.time()-t)
