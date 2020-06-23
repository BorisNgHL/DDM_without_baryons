from matplotlib import pyplot as plt

def read_print(file, comm):
	f = open(file,"r")
	r_arr = []
	DM_arr = []
	
	for line in f:
		arr = []
		arr = f.readline().strip("\n").split("\t")
		if arr != [""]:
			r_arr.append(float(arr[2]))
			DM_arr.append(float(arr[6]))
				
	f.close()
	plt.plot(r_arr,DM_arr,label=comm)
	pass

read_print("T_t=0.txt", "NFW DM density")
read_print("T_result.txt", "Jacky's DM density")
read_print("mod_T_result.txt", "My DM density")
read_print("spline_T_result.txt", "Spline DM density")
plt.xlabel("r (kpc)")
plt.ylabel("DM density (solar mass kpc -3)")
plt.title("DM density against r")
plt.legend()
plt.loglog()
plt.show()