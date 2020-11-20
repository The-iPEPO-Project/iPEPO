
#Physical params
g=2.0  
delta=1


#TEBD params
N_steps_1=10000            
N_steps_2=100000               
N_steps_3=0 #1000000   
N_steps_4=0 #1000000   
N_steps_5=0 #1000000      

dt_1=0.1
dt_2=0.01     
dt_3=0.001   
dt_4=0.0005
dt_5=0.0001

test_interval_1=10     
test_interval_2=50    
test_interval_3=500   
test_interval_4=500   
test_interval_5=500  


#Range of correlations (summary files)
min_sep=0
max_sep=19


#TEBD params
eps=2.0
chi=2        
eta=1.0d-04  


#Environment params
epsENV=1.0d-10   
chiENV=2           
etaENV=1.0d-06   



#Parameter loop vars
xname="g"
x0=${g} 
dx=-0.1
nx=11


#rho_vec initialization params
s1=1.0
s2=1.0



