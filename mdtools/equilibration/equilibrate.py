from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
import mdtraj
from mdtraj.reporters import HDF5Reporter

def equilibrate(topology, positions, temperature=300.0*kelvin):
    """
    Equilibrate an MD system represented by the given topology object
    and provided initial positions in an NPT ensemble

    Parameters
    ----------
    topology : OpenMM.Topology
        Topology object representing the system to simulate
    positions : OpenMM.unit.Quantity([], unit=distance)
        XYZ coordinates of each atom in topology to be used as initial 
        positions
    temperature : OpenMM.unit.Quantity(unit=kelvin)
        Temperature to use for simulation
    """

    # Create System
    ff = ForceField('amber14/protein.ff14SB.xml', 'amber14/tip3p.xml')
    trajtop = mdtraj.Topology.from_openmm(topology)
    system = ff.createSystem(topology, nonbondedMethod=PME, 
                             nonbondedCutoff=1.*nanometer, 
                             constraints=HBonds)

    # Add position restraints that can be tapered off during simulation
    force = CustomExternalForce("k*((x-x0)^2+(y-y0)^2+(z-z0)^2)")
    force.addGlobalParameter("k", 5.0*kilocalories_per_mole/angstroms**2)
    force.addPerParticleParameter("x0")
    force.addPerParticleParameter("y0")
    force.addPerParticleParameter("z0")
    for i in trajtop.select("not water and not (element Na or element Cl) and not element H"):
        force.addParticle(int(i), positions[i].value_in_unit(nanometers))
    system.addForce(force)

    # Setup MD simulation in NPT ensemble
    dt = 0.002*picoseconds
    integrator = LangevinIntegrator(temperature, 1/picosecond, dt)
    barostat   = MonteCarloBarostat(1.0*bar, temperature, 25)
    system.addForce(barostat)
    simulation = Simulation(topology, system, integrator)
    simulation.context.setPositions(positions)

    # Get simulation platform
    platform = Platform.getName(simulation.context.getPlatform())
    print("Running simulation using {}".format(platform))

    simulation.reporters.append(HDF5Reporter('equilibration.h5', 25000))
    simulation.reporters.append(StateDataReporter('equilibration.log', 500, step=True, time=True, volume=True, totalEnergy=True, temperature=True, elapsedTime=True))
    
    # Run minimization
    simulation.minimizeEnergy()

    # Equilibration while tapering off position restraints
    simulation.context.setVelocitiesToTemperature(temperature)
    for i in range(21):
        k = (5.0 - (0.25*i))
        simulation.context.setParameter('k', k*kilocalories_per_mole/angstroms**2)
        simulation.step(25000)

    # Close reporters
    for rep in simulation.reporters:
        rep.close()
        
    return
