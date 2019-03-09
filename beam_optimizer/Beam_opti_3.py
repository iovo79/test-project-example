from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
from math import*
import codecs
from string import Template

""" This script calculates and creates graphs of the deflection in relation
with distance centroid of a grillage of beams for many sections. It also
calculates and graphs the demanded number of beams for the certain spread
length (var:Spread Length) to be covered and the demanded total weight of
the structure in relation with the centroids distance (var:x) in steps
chosen by the user. The length(var:Length) of the beam is also choosable.
The calculations are made for a simple supported beam with a uniform load
(def:simple_beam_full_load).17-10-2014 CMT Copenhagen"""
#-------------------------------------------------
#input
#-------------------------------------------------
#Material
E=210000                            #N/mm2
v=0.3
G=E/(2*(1+v))                       #N/mm2
f_y = 235;f_u = 360                 #N/mm2
#geometry
Spread_length = 8.44                #m Length of area that we have to cover with the beams.
Length = 7.00                       #m
num_of_beams = ()                   #__init__
x=1.00                              #__init__centroid distance between beams starting value
#Create list of sections give name,h,b,tw,tf
Sections =[
#~ ['HEA300',290,300,8.5,14,27],
#~ ['HEA320',310,300,9.0,15.5,27],
#~ ['HEA340',330,300,9.5,16.5,27],
#~ ['HEA360',350,300,10.0,17.5,27],
#~ ['HEA400',390,300,11.0,19.0,27],
#~ ['HEA450',440,300,11.5,21.0,27],
['HEB300',300,300,11.0,19.0,27],
['HEB320',320,300,11.5,20.5,27],
['HEB340',340,300,12.0,21.5,27],
['HEB360',360,300,12.5,22.5,27],
['HEB400',400,300,13.5,24.0,27],
['HEB450',450,300,14.0,26.0,27],
#~ ['IPE300',300,150,7.10,10.7,15],
#~ ['IPE330',330,160,7.50,11.5,18],
#~ ['IPE360',360,170,8.00,12.7,18],
#~ ['IPE400',400,180,8.60,13.5,21],
#~ ['IPE450',450,190,9.40,14.6,21],
#~ ['IPE500',500,200,10.2,16.0,21],\
]
#
Section_props=[]                                                        #__init__ list to save the properties of all section
#Loads definition
#partial uniform load
g1 = 0.00                                                               #__init__ for self weight
q1 = 12.00                                                              #kN/m**2 Scaffolding load
q2 = 1.1*25*0.40                                                        #kN/m**2 Casting concourse slab h=400mm load
load = 0                                                                #__init__ load applied with partial factors
#
out_x = [];out_q = []; out_V = []; out_M = [];                          #__init__lists to be used saving all output
out_w = []; out_1L=[]; out_gw=[];out_n_of_beams=[];I_y =0.00

#FUNCTION DEFINITIONS
def simple_beam_full_load(E=1000*E,L=Length,q=load,I=I_y,x=0.5*Length):
    k=x/L
    V_x = (0.5-k)*q*L
    M_x = (0.5*k*(1-k))*q*L**2
    w_x=1000*k*(1-2*(k**2)+k**3)*((q*L**4)/(24*E*I))
    #~ print "simple beam with uniform load of {0:3.2f}kN/m".format(q)
    #~ print "at point x = {0:3.2f}:\n\
    #~ V(x) = {1:6.3f} kN\n\
    #~ M(x) = {2:6.3f} kN\n\
    #~ δ(x) = {3:6.3f} mm"\
    #~ .format(x,V_x,M_x,w_x)
    #
    R_A = 0.5*q*L
    R_B = R_A
    M_A = 0
    M_B = 0
    w_A = 0
    w_B = 0
    #~ print "Reactions:\n\
    #~ R_A =  {0:6.3f} kN\n\
    #~ R_B =  {1:6.3f} kN\n\
    #~ M_A =  {2:6.3f} kNm\n\
    #~ M_B =  {3:6.3f} kNm\n\
    #~ δ_Α =  {4:6.3f} mm\n\
    #~ δ_B =  {5:6.3f} mm"\
    #~ .format(R_A,R_B,M_A,M_B,w_A,w_B)
    #
    V_max = 0.5*q*L
    x_v_max=0
    M_max=q*(L**2 )/8
    w_max = 1000*5*q*(L**4)/(384*E*I)
    x_M_max = 0.5*L
    #~ print "Max values:\n\
    #~ V_max =  {0:6.3f} kN at x = {1:3.1f}m\n\
    #~ M_max =  {2:6.3f} kNm at x = 0.5L\n\
    #~ δ(max) = {3:6.3f} mm at x = 0.5L\n"\
    #~ .format(V_max,x_v_max,M_max,w_max)
    s_x ={}
    s_x = dict([('V(x)', round(V_x,3)), ('M(x)', round(M_x,3)),\
    ('d(x)', round(w_x,3))])

    S = [round(V_max,3),round(M_max,3),round(w_max,3)]
    return S
#----------------------------------------------------------------

i=0                                                                     #__init__counter
for i in range(len(Sections)):
    #
    #This loop takes the list above of steel sections of type I and
    #returns a list -Section_props[]- which contains lists
    #-Section_props[][]- with the derived properties
    #
    Section = Sections[i][0]
    h = Sections[i][1];
    b = Sections[i][2];
    tw = Sections[i][3];
    tf = Sections[i][4];
    r = Sections[i][5]
    #
    #~ DEBUG BLOCK
    #~ print "Section {0:1.0f}/{1:1.0f} {2:s} {3:3.2f} {4:3.2f} {5:3.2f}\
    #~ {6:3.2f} {7:3.2f}"\
    #~ .format(i,len(Sections),Sections[i][0],Sections[i][1],\
    #~ Sections[i][2],Sections[i][3],Sections[i][4],Sections[i][5])
    #~ DEBUG BLOCK
    #
    d = h-2*tf-2*r                                                        #mm           depth of straight portion of web
    hw = d
    hi = h-2*tf                                                           #mm           inner depth between flanges
    #
    A_c=0.01*(2*tf*b+(h-2*tf)*tw+(4-pi)*(r**2))                           #cm^2         Area section
    I_y = (10**(-8))*0.0001*((1/12)*((b*(h**3))-(b-tw)*((h-2*tf)**3))+\
    0.03*(r**4)+(0.2146*(r**2))*((h-2*tf-0.4468*r)**2))                   #m^4          moment of inertia y-y (strong axis)
    g_w=A_c*(10**(-4))*7850                                               #kgr/m        per meter weight
    #
    #build a list for each section
    #Block for parsing the data
    current_Section = []                                                  #__init__ in every loop
    current_Section.append(Section)
    current_Section.append(I_y)
    current_Section.append(A_c)
    current_Section.append(g_w)
    #
    Section_props.append(current_Section)                                  #save all sections to one list with format [name,I_y,A_c,g_w]


#
#~ DEBUG BLOCK
#~ print len(Section_props)
#~ print Section_props
#~ DEBUG BLOCK
#
print("\n")
j=0                                                                         #__init__counter
for j in range(len(Section_props)):

    print ("------------------------------------------------------------")
    print (" Calculations for Section {0}".format(Section_props[j][0]))
    print ("------------------------------------------------------------")
    x=1.00                                                                  #__reinit__for each section in Section_props list
    current_out_x = []                                                      #__reinit__for each section in Section_props list
    current_out_q = []                                                      #__reinit__for each section in Section_props list
    current_out_V = []                                                      #__reinit__for each section in Section_props list
    current_out_M = []                                                      #__reinit__for each section in Section_props list
    current_out_w = []                                                      #__reinit__for each section in Section_props list
    current_out_1L = []                                                     #__reinit__for each section in Section_props list
    current_out_gw = []                                                     #__reinit__for each section in Section_props list
    current_out_n_of_beams = []                                             #__reinit__for each section in Section_props list


    while x < 2.50:
        x = x +0.10                                                         #loop for every 0.10 m centroid distance
        #
        #~ START DEBUG BLOCK
        print("\n********loop for x={0:3.2f}m*******".format(x))
        g1=0.981*0.01*Section_props[j][3]                                   #kN/m uniform
        load =  1.35*g1 +1.65*q1*x+1.00*q2*x                                #kN/m uniform
        print("DL = {0:3.2f}kN/m | LL1 = {1:3.2f}kN/m**2|LL2 = {2:3.2f}kN/m**2l|Load = {3:3.2f} kN/m"\
        .format(g1,q1,q2,load))
        print ("Section {0} Iy = {1:3.2f}cm^4".format(Section_props[j][0]\
        ,(10**(8))*Section_props[j][1]))
        #~ END DEBUG BLOCK
        #
        s=[]
        s =  simple_beam_full_load(q=load,I=Section_props[j][1])           #Call the function  for calculation of V,M,δ for each case
        #
        #~ START DEBUG BLOCK
        #print output for each x
        print( "x = {0:3.2f}m| Vmax = {1:3.2f}kN| Mmax = {2:3.2f}kNm |δmax = {3:3.2f}mm = L/{4:3.2f}\n"\
        .format(x,s[0],s[1],s[2],round(1000*Length/s[2],0)))
        #~ END DEBUG BLOCK

    #ACTIVE
        #Block for calculating the total number of beams
        #\needed for each x case so total weigth will be derived
        num_of_beams = modf(Spread_length/x)
        if num_of_beams[0] >= 0.5:
            nbeams = ceil(Spread_length/x)
        else:
            nbeams = floor(Spread_length/x)
        print( "{0:3.0f} of {1} of total weight of {2:3.2f} Kgr"\
        .format(nbeams,Section_props[j][0],Section_props[j][3]*nbeams*Length))
    #ACTIVE

    #
    #build a list for each x
    #Block for parsing the data
        current_out_x.append(x)
        current_out_q.append(load)
        current_out_V.append(s[0])
        current_out_M.append(s[1])
        current_out_w.append(s[2])
        current_out_1L.append(1000*Length/s[2])
        current_out_n_of_beams.append(nbeams)
        current_out_gw.append(Section_props[j][3]*nbeams*Length)
    #
    #save all output for each x to one list for each Section in a form of output for each section
    #for example if you have 7 sections and 10 steps for each x calculation  the list
    # out_w (max deformation) would be [ zero[w0,...,w9]....6th[w0,...,w9]]
    out_x.append(current_out_x)                                         #xi centroid distance
    out_q.append(current_out_q)                                         #kN/m load applied (total)
    out_V.append(current_out_V)                                         #kN  Shear max value
    out_M.append(current_out_M)                                         #kNm Moment max value
    out_w.append(current_out_w)                                         #mm deformation max value
    out_1L.append(current_out_1L)                                       #deformation max value in terms of Length
    out_n_of_beams.append(current_out_n_of_beams)                       #number of beams needed for the current (in loop) x[i] centroid value
    out_gw.append(current_out_gw)                                       #total steel weight deployed for the current (in loop) x[i] centroid value
#
#~ -START--DEBUG BLOCK
#~ to check the output lists.
#~ k = 0                                                                #__init__counter
#~ for k in range(len(out_x)):
    #~ print "loop {0:1.0f}/{0:1.0f}".format(k,len(out_x))
    #~ print "Section {0}".format(Section_props[k][0])
    #~ print out_x[k]
    #~ print out_1L[k]
#~ -END DEBUG BLOCK

#PLOTTING
k = 0
plt.figure(1)
plt.suptitle('SOLUTION WITH {0} PROFILE'\
.format(Section_props[k][0][0:3]), fontsize=20)

plt.subplot(231)
for k in range(len(out_x)):
    plt.plot(out_x[k],out_1L[k],lw=2,label=Section_props[k][0])
#plt.xlim([0.00,2.6])
plt.xlabel('$d_{i}$',fontsize=20)
plt.ylabel('L/$\delta$',fontsize=20)
#plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,ncol=3,\
#mode="expand", borderaxespad=0.,fontsize='xx-small')
plt.legend(fontsize='large')
plt.grid(True)

plt.subplot(234)
for k in range(len(out_x)):
    plt.plot(out_x[k],out_w[k],lw=2,)
plt.xlabel('$d_{i}$',fontsize=20)
plt.ylabel('$\delta$ [mm]',fontsize=20)
plt.legend(fontsize='large')
plt.grid(True)

plt.subplot(232)
for k in range(len(out_x)):
    plt.plot(out_x[k],out_M[k],lw=2)
plt.xlabel('$d_{i}$',fontsize=20)
plt.ylabel('$M_{max}$[kNm]',fontsize=20)
plt.grid(True)

plt.subplot(235)
for k in range(len(out_x)):
    plt.plot(out_x[k],out_V[k],lw=2)
plt.xlabel('$d_{i}$',fontsize=20)
plt.ylabel('$V_{max}$ [kN]',fontsize=20)
plt.grid(True)

plt.subplot(233)
for k in range(len(out_x)):
    plt.plot(out_x[k],out_n_of_beams[k],lw=2)
plt.xlabel('$d_{i}$',fontsize=20)
plt.ylabel('No of Beams',fontsize=20)
plt.grid(True)

plt.subplot(236)
for k in range(len(out_x)):
    plt.plot(out_x[k],out_gw[k],lw=2)
plt.xlabel('$d_{i}$',fontsize=20)
plt.ylabel('[kgr]',fontsize=20)
plt.grid(True)

plt.show()
