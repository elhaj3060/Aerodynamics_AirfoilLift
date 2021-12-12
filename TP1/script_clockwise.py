counterclockwise = open('NACA1410.dat' , 'r')
clockwise = open("CW_NACA1410.dat", "w")
a = counterclockwise.readlines()


for i in range(len(a)-1, -1,-1):
	print(i)
	print(a[i])
	clockwise.write(a[i])

counterclockwise.close()
clockwise.close()
 
