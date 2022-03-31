#!/usr/local/bin/python

#############################################################
# Simple MD of Lennard Jones charged or uncharged particles #
# Alexandre Bonvin, Aalt Jan van Dijk, Utrecht University   #
# Updated to a singular script with a modern look by Douwe  #
# Schulte, Utrecht University (2022)                        #
#                                                           #
# Adapted from a script from Patrick Fuchs, Uni. Paris VI   #
#############################################################


##################
# Import modules #
##################
from math import sqrt,log,sin,cos
from random import random,seed
from enum import Enum
from tkinter import Tk, Canvas, DoubleVar, StringVar
from tkinter.ttk import Label, Button, Style, Frame, Notebook, Entry
import sys


#####################
# Define parameters #
#####################

nAtoms       = 20                    # Number of atoms
Radius       = 25.0                  # Beware that Radius must be in a good range (according to nAtoms)
                                     # In order to be able to place all atoms
Mass         = 10.0                  # Atom mass
Rmin         = 2.24 * Radius         # Distance at which Rmin is minimal
BoxDim       = [500,500]             # Box dimension
Atom_Coord   = []                    # List of the form : [nAtoms][2]
Epsilon      = 2 * Radius            # Well depth
Dielec       = 1.0                   # Dielectric constant
qat          = 2 * Radius            # Atom absolute charge
frac_neg     = 0.5                   # Fraction negative charges
OverlapFr    = 0.0                   # Fraction of overlap allowed
CutOff       = 250                   # Non-bonded cutoff
CutOffSquare = CutOff**2             # Precalculated square
speed        = 5                     # Canvas update speed
cstboltz     = 0.00198722            # Boltzmann's constant in kcal/mol/K
cstboltz     = 1000*cstboltz/4.18    # In kJ/mol/K
Seed         = 42                    # Random number seed
# Steepest Descent parameters
drinit       = 1.00	                 # dr from EM
drmin        = 0.00001               # Minimum dr value to step EM
drmax        = 5.00                  # Maximum dr
alpha        = 1.05                  # Scaling factor for dr if Enew < Eold
beta         = 0.90                  # Scaling factor for dr if Enew > Eold
deltaE       = 0.001                 # Energy difference threshold to stop EM
normFmin     = 0.001                 # Minimum force norm to step EM
# Verlet parameters
Temperature  = 300.0                 # Temperature in K
timestep     = 5.0E-3	             # MD time step
# Set specific behaviour for practical session 4
CalculateEnergyPeriodic = True       # Practical #4 part 1
ShowOtherEnergyCutoffResults = False # Practical #4 part 2

# Additional program specific parameters
Minimizers = Enum("Minimisers", "SteepestDescent Verlet")
Minimizer = Minimizers.Verlet
drstep = drinit
Iterations = 0
canvas_event = None
Color = []
ATOM = []


##############################
# Steepest descent minimizer #
##############################

def steepest_descent(atom_coord,drstep,forces):
    """
    This function gets as input parameters:
    - atom_coord, a vector containing the x and y position and the charge of the i atoms
    - drstep, the displacement for the minimizer
    - force, a vector containing the x and y components of the force on the atoms

    The function returns a list array (vector containing the new positions)
    Implement in the following loop over all atoms the using the steepest descent algorithm
    A few hints:
    - powers in python are given by **, e.g.: x to the square is x**2
    - squared root x: sqrt(x)
    - avoid dividing by zero
    """
    new_positions=[]

    # 1) First calculate the norm of the total force vector
    normf = 0.0
    for force in forces:
        normf=normf+force[0]**2.0+force[1]**2.0
    normf=sqrt(normf)

    if normf < 0: return atom_coord, normf

    # 2) Then move the particles
    for (coord, force) in zip(atom_coord, forces):
        r0x=coord[0]		# Coordinates
        r0y=coord[1]

# Insert below the lines defining the new coordinates based on the old ones + forces + drstep.
#
# Forces are contained in force[0] for the x force component and force[1] for the y force.
# The step size for the move is given by drstep.
#
# ====>>>>>
        sx=force[0]/normf
        sy=force[1]/normf
        r0xnew=r0x+drstep*sx
        r0ynew=r0y+drstep*sy
# <<<<<====

        new_positions.append([r0xnew,r0ynew,coord[2]])
    return new_positions,normf


#####################
# Verlet integrator #
#####################

def verlet(atom_coord,forces,dtstep,old_atom_coord,mass):
    """
    This function gets as input parameters:
    - `atom_coord`, a vector containing the x and y position and the charge of the i atoms
    - `old_atom_coord`, a vector containing the x and y positions from the previous MD step
    - `forces`, a vector containing the x and y components of the force on the atoms
    - `dtstep`, the integration time step
    The function returns a list containing the new positions.

    Implement in the following loop between the arrows the Verlet MD algorithm.
    A few hints:
    - Powers in python are given by **, e.g.: x to the square is `x**2`
    - Squared root x: `sqrt(x)`
    - Indents are important in python
    """

    new_positions=[]
    for coord,old_coord,force  in zip(atom_coord, old_atom_coord, forces):
        r0x=coord[0]			# Coordinates
        r0y=coord[1]
        old_r0x=old_coord[0]		# Old coordinates
        old_r0y=old_coord[1]

# Insert below the lines defining the new x and y positions based on the old ones + forces + mass + dtstep.
#
# Forces are contained in force[0] for the x force component and force[1] for the y force.
# The step size for the move is given by dtstep.
#
# ====>>>>>
        new_r0x = ...
        new_r0y = ...
# <<<<<====
        new_positions.append([new_r0x,new_r0y,coord[2]])
    return new_positions


def verlet_step1(atom_coord,velocity,forces,dtstep,mass):
    """The first step for Verlet"""
    global Ene,EneLJ,EneCoul,ELJ2,ELJ4
    new_positions=[]
    for coord, vel, force in zip(atom_coord, velocity, forces):
        r0x=coord[0]+dtstep*vel[0]+0.5*dtstep**2*force[0]/mass
        r0y=coord[1]+dtstep*vel[1]+0.5*dtstep**2*force[1]/mass
        new_positions.append([r0x,r0y,coord[2]])
    Ene,EneLJ,EneCoul,ELJ2,ELJ4 = calculate_energy(Atom_Coord,Epsilon,Rmin,Dielec,CutOffSquare,BoxDim)
    return new_positions

def calculate_velocities(old_atom_coord,atom_coord,dtstep):
    """Calculate velocities based on old and new positions"""
    velocities=[]
    for coord, old_coord in zip(atom_coord, old_atom_coord):
        v0x=(coord[0]-old_coord[0])/(2*dtstep)
        v0y=(coord[1]-old_coord[1])/(2*dtstep)
        velocities.append([v0x,v0y])
    return velocities


##########################
# Move particles with MD #
##########################

def simulate():
    """Execute the simulation"""
    global Atom_Coord,Radius,Mass,BoxDim,Epsilon,Rmin,CutOffSquare,Iterations,Ene,Old_Atom_Coord
    global Velocity,timestep,report_var_total,report_var_subenergies, drstep, Ene_prev
    global Color,report_var_time,Dielec,root,atom_canvas,speed,canvas_event

    Force = calculate_force(Atom_Coord,Epsilon,Rmin,Dielec,CutOffSquare,BoxDim)
    tmp=Atom_Coord

    if Minimizer == Minimizers.SteepestDescent:
        if Iterations == 0: Ene_prev=Ene
        Atom_Coord, normF=steepest_descent(Atom_Coord,drstep,Force)

    if Minimizer == Minimizers.Verlet:
        if Iterations == 0:
            Old_Atom_Coord=Atom_Coord
            Atom_Coord=verlet_step1(Atom_Coord,Velocity,Force,timestep,Mass)

        Atom_Coord=verlet(Atom_Coord,Force,timestep,Old_Atom_Coord,Mass)
        Velocity=calculate_velocities(Old_Atom_Coord,Atom_Coord,timestep)

    Old_Atom_Coord=tmp
    Ene,EneLJ,EneCoul,ELJ2,ELJ4 = calculate_energy(Atom_Coord,Epsilon,Rmin,Dielec,CutOffSquare,BoxDim)
    Kin,temperature=calculate_temperature(Velocity,nAtoms,cstboltz,Mass)

    # Update drstep
    if Minimizer == Minimizers.SteepestDescent:
        if Ene < Ene_prev:
            drstep = min(drmax, drstep * alpha)
        else:
            drstep = drstep * beta

    # Update top labels
    report_var_time.set("Step: %d Time: %8.3f" % (Iterations,float(Iterations)*timestep))
    report_var_total.set("Etot: %6.1f Ekin: %6.1f Epot: %6.1f" % (Ene+Kin,Kin,Ene))
    if ShowOtherEnergyCutoffResults:
        report_var_subenergies.set("Elj: %6.2f Elj2: %6.2f Elj4: %6.2f Ecoul: %6.1f Temp: %6.1f" % (EneLJ,ELJ2,ELJ4,EneCoul,temperature))
    else:
        report_var_subenergies.set("Elj: %6.1f Ecoul: %6.1f Temp: %6.1f" % (EneLJ,EneCoul,temperature))

    # Apply boundary conditions
    for coord, old_coord in zip(Atom_Coord, Old_Atom_Coord):
        for i in range(2): # i=0 -> case x coordinate ; i=1 -> case y coordinate
            if coord[i] < 0:
                coord[i] += BoxDim[i]
                old_coord[i] += BoxDim[i]
            if coord[i] > BoxDim[i]:
                coord[i] -= BoxDim[i]
                old_coord[i] -= BoxDim[i]

    # Draw new canvas coordinates
    for atom, coord in zip(ATOM, Atom_Coord):
        x, y = coord[0], coord[1]
        atom_canvas.coords(atom, x + Radius, y + Radius, x - Radius, y - Radius)

    # Print to terminal window
    if Iterations % 20 == 0:
        if ShowOtherEnergyCutoffResults:
            print("Step: %4d Time: %8.3f Etot: %6.1f Ekin: %6.1f Epot: %6.1f Elj: %6.1f Elj2: %6.1f Elj4: %6.1f Ecoul: %6.1f Temp: %6.1f" % (Iterations,float(Iterations)*timestep,Ene+Kin,Kin,Ene,EneLJ,ELJ2,ELJ4,EneCoul,temperature))
        else:
            print("Step: %4d Time: %8.3f Etot: %6.1f Ekin: %6.1f Epot: %6.1f Elj: %6.1f Ecoul: %6.1f Temp: %6.1f" % (Iterations,float(Iterations)*timestep,Ene+Kin,Kin,Ene,EneLJ,EneCoul,temperature))

    # Stopping conditions
    if Minimizer == Minimizers.SteepestDescent and (abs(Ene - Ene_prev) < deltaE or drstep < drmin or normF < normFmin):
        print("STOPPING... deltaE <",deltaE,", or drstep <",drmin,", or normF <",normFmin)
        outtext="Step: %4d Epot: %6.1f Elj: %6.1f Ecoul: %6.1f deltaE: %10.6f <normF>: %8.6f dr: %8.6f" % (Iterations,Ene,EneLJ,EneCoul,Ene - Ene_prev,normF,drstep)
        print(outtext)
    elif temperature > 1000000:
        print("The system is exploding !!!")
        print("Step: %4d Time: %8.3f" % (Iterations,float(Iterations)*timestep))
        print("Etot: %6.1f Ekin: %6.1f Epot: %6.1f" % (Ene+Kin,Kin,Ene))
        print("Elj: %6.1f Ecoul: %6.1f Temp: %6.1f" % (EneLJ,EneCoul,temperature))
        print("Emergency stop")
    else:
        Ene_prev=Ene
        Iterations=Iterations+1
        canvas_event=atom_canvas.after(speed,simulate)


####################
# Energy functions #
####################

# Calculate Lennard Jones from the squared distance
def LJ2(r2, epsilon, sigma6):
    # Uncomment the following lines to get a more obvious mathematical implementation
    # r = sqrt(r2)
    # sigma = sigma6**(1/6)
    # return epsilon*((sigma/r)**12 - (sigma/r)**6)

    # The following implementation is significantly faster so this is the default
    Z = (1/r2)**3 * sigma6
    return epsilon * Z * (Z-1)

# Classical Coulomb from the squared distance
def Coulomb2(r,dielec,qa,qb):
    return qa*qb/(dielec*sqrt(r))

# Calculate energy Evdw + Ecoulomb (used squared distance)
def calculate_energy(coord,epsilon,rmin,dielec,cutoffsquare,boxdim,elec=1):
    global CalculateEnergyPeriodic
    cutoff2=2.0*rmin; cutoff2sq=cutoff2**2
    cutoff4=4.0*rmin; cutoff4sq=cutoff4**2
    Ene = 0.0; distsquare = 0
    ELJ = 0.0; ECoul=0.0
    ELJ2 = 0.0; ELJ4 = 0.0
    rmin_exp6 = rmin**6
    # Doubly nested loop over all particle pairs
    for i in range(len(coord)-1):
        for j in range(i+1,len(coord)):
            # Calculate the squared atomic distance
            distsquare = 0
            for k in range(2):
                tmp = coord[j][k] - coord[i][k]
                # Chooses the nearest image
                if CalculateEnergyPeriodic:
                    halfbox = boxdim[k]/2
                    tmp = tmp - SignR(halfbox,tmp-halfbox) - SignR(halfbox,tmp+halfbox)
                distsquare += tmp**2
            # Compute vdw and Coulomb energy
            if distsquare < cutoffsquare:
                qa = coord[i][2]
                qb = coord[j][2]
                vdw  = LJ2(distsquare, epsilon, rmin_exp6)
                Ene += vdw
                ELJ += vdw
                if elec:
                    CC = Coulomb2(distsquare,dielec,qa,qb)
                    Ene+=CC
                    ECoul+=CC
                if distsquare < cutoff4sq:
                    ELJ4 += vdw
                    if distsquare < cutoff2sq:
                        ELJ2 += vdw

    return Ene,ELJ,ECoul,ELJ2,ELJ4

# Calculate kinetic energy and temperature
def calculate_temperature(vel,nat,k,mass):
    v2=0.0
    for velocity in vel:
        v2=v2+velocity[0]**2+velocity[1]**2
    nkt=0.5*mass*v2		# Kinetic energy equals 0.5*m*v**2
    kin=v2*0.5*mass
    temp=nkt/(nat*k)		# N*k*T=Kinetic Energy
    return kin,temp


###################
# Force functions #
###################

# Force LJ (use squared distance)
def calculate_lennard_jones(distsquare, epsilon, rmin_exp6,xi):
    rij=sqrt(distsquare)
    Z = (1/distsquare)**3 * rmin_exp6
    dedz=epsilon*(2*Z-1)
    dzdr=rmin_exp6*(-6.0/rij**(7.0))
    drdx=xi/rij
    return dedz*dzdr*drdx

# Force Coulomb (use squared distance)
def calculate_coulomb(distsquare,dielec,qa,qb,xi):
    rij=sqrt(distsquare)
    dedr=-1.0*(qa*qb/dielec)*(1/distsquare)
    drdx=xi/rij
    return dedr*drdx

# Calculate force from Evdw + Ecoulomb (uses squared distance)
def calculate_force(coord,epsilon,rmin,dielec,cutoffsquare,boxdim):
    Force=[] ; distsquare = 0
    rmin_exp6 = rmin**6
    # Doubly nested loop over all particle pairs
    for i in range(len(coord)):
        tmpforce=[0.0,0.0]
        for j in range(len(coord)):
            if not i==j:
                # Calculate the squared atomic distance
                distsquare = 0
                for k in range(2):
                    tmp = coord[j][k] - coord[i][k]
                    # Chooses the nearest image
                    halfbox = boxdim[k]/2
                    tmp = tmp - SignR(halfbox,tmp-halfbox) - SignR(halfbox,tmp+halfbox)
                    distsquare += tmp**2
                # Compute vdw force
                if distsquare < cutoffsquare:
                    qa = coord[i][2]
                    qb = coord[j][2]
                    for k in range(2):
                        tmp = coord[j][k] - coord[i][k]
                        ff = calculate_lennard_jones(distsquare, epsilon, rmin_exp6,tmp)
                        ff += calculate_coulomb(distsquare,dielec,qa,qb,tmp)
                        tmpforce[k]+=ff
        Force.append(tmpforce)
    return Force


###################
# Other functions #
###################

# Normal Distance
def dist(A,B):
    return sqrt((A[0]-B[0])**2+(A[1]-B[1])**2)

# Change sign
def SignR(a,b):
    if b > 0:
        return a
    else:
        return -a

# Color particules based on charge
def charge_color(charge,qat):
    if charge == qat:
        return "white"
    else:
        return "#333333"


##################
# Initialization #
##################

# Generate random coordinates
def InitConf(n,dim,radius,qat,frac_neg):
    global Seed
    seed(Seed)
    print("Initializing box, please wait...", end='')
    tmp_coord = []
    ntrial = 0
    i = 1
    # Fix first atom
    x = random()*(dim[0]-2*radius)+radius
    y = random()*(dim[1]-2*radius)+radius
    nneg = int(float(n) * frac_neg)
    charge = -qat
    if nneg == 0: charge = qat
    tmp_coord.append([x,y,charge])
    for negative in [-1, 1]:
        while negative == -1 and i < nneg or negative == 1 and i < n:
            x = random()*(dim[0]-2*radius)+radius
            y = random()*(dim[1]-2*radius)+radius
            # Check wether the new particle overlaps with an existing one
            OVERLAP = False
            for j in range(i):
                if dist(tmp_coord[j],[x,y]) < (1-OverlapFr)*2*radius:
                    OVERLAP = True
            if not OVERLAP:
                tmp_coord.append([x,y,negative * qat])
                i += 1
            ntrial = ntrial + 1
            if ntrial > 100000:
                print('error')
                print("Initialisation failed")
                print("==> Reduce radius or number of atoms")
                sys.exit()
    print("done")
    return tmp_coord

# Generate random charges
def InitCharge(n,qat,frac_neg):
    global Atom_Coord
    print("Initializing charges, please wait...", end='')
    i = 0
    nneg = int(float(n) * frac_neg)
    charge = -qat
    if nneg == 0: charge = qat
    Atom_Coord[i][2]=charge
    i += 1
    while i < nneg:
        Atom_Coord[i][2]=-qat
        i += 1
    while i < n:
        Atom_Coord[i][2]=qat
        i += 1
    print("done")

# Generates initial velocities according to Maxwell distribution
def InitVel(n,temperature,cstboltz,mass):
    global Seed
    seed(Seed)
    stdev=sqrt(cstboltz*temperature/mass)
    print("Initializing velocities, please wait...", end='')
    tmp_vel=[]
    for i in range(n):
        # Generate random numbers according to Gaussian:
        r1=random()
        r2=random()
        x1=sqrt(-2.0*log(r1))*cos(r2)*stdev
        x2=sqrt(-2.0*log(r1))*sin(0.5*r2)*stdev
        tmp_vel.append([x1,x2])
    # Remove overall motion
    vxt=0.0
    vyt=0.0
    for item in tmp_vel:
        vxt+=item[0]
        vyt+=item[1]
    for item in tmp_vel:
        item[0] -= vxt/float(n)
        item[1] -= vyt/float(n)
    # Scaling factor is used to get temperature exactly equal to desired temperature
    kin,tt=calculate_temperature(tmp_vel,n,cstboltz,mass)
    scaling=sqrt(temperature/tt)
    vel=[]
    for item in tmp_vel:
        vx=item[0]*scaling
        vy=item[1]*scaling
        vel.append([vx,vy])
    print("done")
    return vel

########################################
# Various functions for input + layout #
########################################

# Setup system
def set_up_atoms(repack=1):
    global Iterations,Velocity,Temperature,Mass,cstboltz,atom_canvas,ATOM,Atom_Coord,Color
    ATOM = []

    if repack==1:
        Atom_Coord = InitConf(nAtoms,BoxDim,Radius,qat,frac_neg)
        Color = []
        for i in range(nAtoms):
            Color.append(charge_color(Atom_Coord[i][2],qat))
        Velocity=InitVel(nAtoms,Temperature,cstboltz,Mass)

    if repack==2:
        InitCharge(nAtoms,qat,frac_neg)
        Color = []
        for i in range(nAtoms):
            Color.append(charge_color(Atom_Coord[i][2],qat))

    for (color, atom) in zip(Color, Atom_Coord):
        x, y = atom[0], atom[1]
        ATOM.append(atom_canvas.create_oval(x + Radius,y + Radius,x - Radius,y - Radius,fill=color))
    update_energy()

# Set number of particles
def set_r(event):
    global nAtoms
    nAtoms=int(r.get())
    update_canvas()

# Set atom Radius
def set_size(event):
    global Radius,Rmin
    Radius=int(size.get())
    Rmin = 2 * Radius
    update_canvas()

# Set epsilon for Lennard-Jones
def set_vdw1(event):
    global Epsilon
    Epsilon=int(vdw1.get())
    update_canvas(0)

# Set sigma for Lennard-Jones
def set_vdw2(event):
    global Rmin
    Rmin=int(vdw2.get())
    update_canvas(0)

# Set charge fraction
def set_frac(event):
    global frac_neg
    frac_neg=float(frac.get())
    update_canvas(2)

# Set particle charge
def set_q(event):
    global qat
    qat=float(q.get())
    update_canvas(2)

# Set dielectric constant
def set_diel(event):
    global Dielec
    Dielec=float(diel.get())
    update_canvas(0)

# Set Temperature
def set_temp(event):
    global Temperature
    Temperature=float(temp.get())
    update_canvas(0)

def set_tstep(event):
    global timestep,Velocity,nAtoms,Temperature,cstboltz,Mass
    timestep=float(tstep.get())
    update_canvas(0)
    Velocity=InitVel(nAtoms,Temperature,cstboltz,Mass)

# Set minimum Force norm difference for stop condition
def set_dFmin(event):
    global normFmin
    normFmin=float(Fmin.get())
    update_canvas(0)

# Set minimum Energy difference for stop condition
def set_deltaE(event):
    global deltaE
    deltaE=float(Emin.get())
    update_canvas(0)

# Set initial displacement for minimizer
def set_dxstep(event):
    global drinit
    drinit=float(dxstep.get())
    update_canvas(0)

# Set alpha factor for increasing dr
def set_alpha(event):
    global alpha
    alpha=float(alphafactor.get())
    update_canvas(0)

# Set beta factor for decreasing dr
def set_beta(event):
    global beta
    beta=float(betafactor.get())
    update_canvas(0)

# Update energy
def update_energy():
    global Atom_Coord,BoxDim,Epsilon,Rmin,CutOffSquare,Iterations,Ene
    global Dielec,Velocity,Mass,cstboltz

    Ene,EneLJ,EneCoul,ELJ2,ELJ4 = calculate_energy(Atom_Coord,Epsilon,Rmin,Dielec,CutOffSquare,BoxDim)
    Kin,temperature=calculate_temperature(Velocity,nAtoms,cstboltz,Mass)

    report_var_time.set("Step: %d Time: %8.3f" % (Iterations,float(Iterations)*timestep))
    report_var_total.set("Etot: %6.1f Ekin: %6.1f Epot: %6.1f" % (Ene+Kin,Kin,Ene))
    if ShowOtherEnergyCutoffResults:
        report_var_subenergies.set("Elj: %6.2f Elj2: %6.2f Elj4: %6.2f Ecoul: %6.1f Temp: %6.1f" % (EneLJ,ELJ2,ELJ4,EneCoul,temperature))
    else:
        report_var_subenergies.set("Elj: %6.1f Ecoul: %6.1f Temp: %6.1f" % (EneLJ,EneCoul,temperature))


def update_canvas(repack=1):
    global Iterations, atom_canvas, Atom_Coord, ATOM, Color
    atom_canvas.delete("all")
    set_up_atoms(repack)
    update_energy()
    Iterations = 0

################
# MAIN PROGRAM #
################

def die():
    sys.exit()

def select_minimizer(*args):
    global Minimizer, minimizer_selector
    reset()
    if minimizer_selector.index('current') == 0: # First tab is Steepest Descent
        Minimizer = Minimizers.SteepestDescent
        selected_method_text.set("Active method: Steepest Descent")
    else:
        Minimizer = Minimizers.Verlet
        selected_method_text.set("Active method: Verlet")

def start():
    global canvas_event
    if canvas_event==None:
        simulate()

def stop():
    global canvas_event
    if canvas_event != None: atom_canvas.after_cancel(canvas_event)
    canvas_event = None
    update_canvas(0)

def reset():
    global canvas_event
    if canvas_event != None: atom_canvas.after_cancel(canvas_event)
    canvas_event = None
    update_canvas()

root = Tk()
root.winfo_toplevel().title("MolMod Practical")
root.bind("<Escape>", die)
root.bind('<Control-c>', die)

top=Frame(root)
top.pack(side='top')
title=Frame(top)
title.pack(side='top')
labels=Frame(top)
labels.pack(side='top')
buttons=Frame(top)
buttons.pack(side='bottom')
atom_canvas=Canvas(root, width=BoxDim[0], height=BoxDim[1],bg="#ccddff")
atom_canvas.pack()
minimizer_selector=Notebook(root)
minimizer_selector.pack(side='bottom')
steepest_descent_pack=Frame(minimizer_selector)
steepest_descent_pack.pack(side='top')
verlet_pack=Frame(minimizer_selector)
verlet_pack.pack(side='top')
Style().configure("Notebook", foreground="black")
minimizer_selector.add(steepest_descent_pack, text="Steepest Descent")
minimizer_selector.add(verlet_pack, text="Verlet")
minimizer_selector.bind("<<NotebookTabChanged>>", select_minimizer)
minimizer_selector.select(1)
selected_method=Frame(root)
selected_method.pack(side='bottom')
low2=Frame(root)
low2.pack(side='bottom')
low1=Frame(root)
low1.pack(side='bottom')

r=DoubleVar()
size=DoubleVar()
vdw1=DoubleVar()
vdw2=DoubleVar()
frac=DoubleVar()
diel=DoubleVar()
q=DoubleVar()
temp=DoubleVar()
tstep=DoubleVar()
Emin=DoubleVar()
Fmin=DoubleVar()
alphafactor=DoubleVar()
betafactor=DoubleVar()
q=DoubleVar()
temp=DoubleVar()
dxstep=DoubleVar()

# Create an entry with a label
def create_entry(pack, text, var, bound_var, callback):
    Label(pack,text=text + " =").pack(side='left')
    var.set(bound_var)
    tstep_entry=Entry(pack,width=6,textvariable=var)
    tstep_entry.pack(side='left')
    tstep_entry.bind('<Return>', callback)
    tstep_entry.bind('<FocusOut>', callback)

# Set up the general parameters
create_entry(low1, "Atoms", r, nAtoms, set_r)
create_entry(low1, "VDW radius", size, Radius, set_size)
create_entry(low1, "VDW ε", vdw1, Epsilon, set_vdw1)
create_entry(low1, "VDW σ", vdw2, Rmin, set_vdw2)
create_entry(low2, "Coulomb param: fraction negative", frac, frac_neg, set_frac)
create_entry(low2, "Charge", q, qat, set_q)
create_entry(low2, "Dielec", diel, Dielec, set_diel)

# Steepest Descent Paramaters
create_entry(steepest_descent_pack, "DeltaE threshold", Emin, deltaE, set_deltaE)
create_entry(steepest_descent_pack, "dFmin", Fmin, normFmin, set_dFmin)
create_entry(steepest_descent_pack, "dr init", dxstep, drinit, set_dxstep)
create_entry(steepest_descent_pack, "α", alphafactor, alpha, set_alpha)
create_entry(steepest_descent_pack, "β", betafactor, beta, set_beta)

# Verlet Parameters
create_entry(verlet_pack, "T (K)", temp, Temperature, set_temp)
create_entry(verlet_pack, "Timestep", tstep, timestep, set_tstep)

# Set up title
Label(title,text="EM & MD",foreground='blue',font='times 18 bold').pack(side='left')

# Set up reporting labels
report_var_time = StringVar()
Label(labels,textvariable=report_var_time).pack(side='top')

report_var_total = StringVar()
Label(labels,textvariable=report_var_total).pack(side='top')

report_var_subenergies = StringVar()
Label(labels,textvariable=report_var_subenergies).pack(side='top')

selected_method_text = StringVar()
Label(selected_method,textvariable=selected_method_text).pack(side='top')

# Set up buttons
Style().configure("TButton", padding=1, relief="flat")
Style().configure("Start.TButton", foreground='blue')
Style().configure("Stop.TButton", foreground='red')
Style().configure("Reset.TButton", foreground='green')

Button(buttons,text='Start',command=start,style="Start.TButton").pack(side='left',fill='x')
Button(buttons,text='Stop',command=stop,style="Stop.TButton").pack(side='left')
Button(buttons,text='Reset',command=reset,style="Reset.TButton").pack(side='left')

# Set up the positions of the atoms and start the simulation
set_up_atoms()

print("Click on 'Start' to go ahead")
print("Use <ESC> or 'X' to quit")

root.mainloop()
