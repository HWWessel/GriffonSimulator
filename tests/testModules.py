# Test the different modules separately with a basic test case as part of the integration of the modules.
# Author/s: Jose Felix Zapata Usandivaras
# Date: 30/12/2018

# --------------------------- IMPORT MODULES ---------------------------

from DataLayer.JsonInterpreter import JsonInterpreter           # Import the json interpreter
import CombustionModule.Geometries as Geom                      # Import the Geometry module
import CombustionModule.Nozzle as Noz                           # Import the Nozzle module
from CombustionModule.Combustion import CombustionObject        # Import the CombustionObject
from MassEstimationModule.system import System                  # Import the system class
from TrajectoryModule.Drag import *                             # Import the Drag library
from TrajectoryModule.Density import DensityLaw                 # Import the density-law library
from TrajectoryModule.Trajectory import TrajectoryObject        # Import the trajectory object
import numpy as np                                              # Import numpy

# ------------------------ FUNCTION DEFINITIONS ------------------------


def generate_data_layer(data_file="Thermodynamic Data 36 bar OF 0,1 to 8,0 H2O2 87,5.json"):
    """ generate_data_layer instantiates all of the objects that
    are required as data layer for the program to run properly.
    :param data_file: name of the json file to be used as reference for the simulation (str)
    :return JsonInterpreter """

    # Pass the file name
    data_directory = "C:/Users/Felix Zapata/Google Drive/" \
                     "MSc in AE ISAE - SUPAERO/Griffon Project/GriffonSimulator/data"

    # Return the output
    return JsonInterpreter(file_name="/".join([data_directory, data_file]))


def test_combustion():
    """ perform the test over the combustion module """

    # ------------ Define parameters:

    geometric_params = {'L': 0.5, 'rintInitial': 0.04, 'rext0': 0.1}

    nozzle_params = {'At': 0.000589, 'expansion': 5.7, 'lambda_e': 0.98, 'erosion': 0}

    design_params = {'gamma': 1.27, 'p_chamber': 4000000, 'p_exit': 100000,
                     'c_star': 1500, 'isp': 230, 'thrust': 30000}

    simulation_params = {'ox_flow': 6, 'safety_thickness': 0.005, 'dt': 0.05}

    # ------------- Generate objects:

    geometry_obj = Geom.OneCircularPort(**geometric_params)
    nozzle_obj = Noz.Nozzle(**nozzle_params)
    nozzle_obj.set_design(**design_params)
    json_interpreter = generate_data_layer()

    # Instantiate the combustion module
    combustion_obj = CombustionObject(json_interpreter=json_interpreter,
                                      geometry_object=geometry_obj,
                                      nozzle_object=nozzle_obj)

    # -------------- Run simulation & Plot:

    combustion_obj.run_simulation_constant_fuel_sliver(**simulation_params)

    # Print the module
    print(combustion_obj)

    # Plot the results
    combustion_obj.plot_results()


def test_combustion_onera_data():
    """ perform the test over the combustion module but with Onera Test Data """

    # ------------ Define parameters:

    geometric_params = {'L': 0.235, 'rintInitial': 0.0186658954, 'rext0': 0.2}

    nozzle_params = {'At': 0.000038, 'expansion': 6.3, 'lambda_e': 0.98, 'erosion': 0}

    design_params = {'gamma': 1.27, 'p_chamber': 4000000, 'p_exit': 100000,
                     'c_star': 1500, 'isp': 230, 'thrust': 30000}

    simulation_params = {'ox_flow': 0.0855, 'safety_thickness': 0.005, 'dt': 0.05, 'max_burn_time': 8.25}

    # ------------- Generate objects:

    geometry_obj = Geom.OneCircularPort(**geometric_params)
    nozzle_obj = Noz.Nozzle(**nozzle_params)
    # nozzle_obj.set_design(**design_params)                            # Do not set the design of the nozzle
    json_interpreter = generate_data_layer(data_file="Thermodynamic Data Onera 41 bar H2O2 87_5.json")

    # Instantiate the combustion module
    combustion_obj = CombustionObject(json_interpreter=json_interpreter,
                                      geometry_object=geometry_obj,
                                      nozzle_object=nozzle_obj)

    # -------------- Run simulation & Plot:

    combustion_obj.run_simulation_constant_fuel_sliver(**simulation_params)

    # Print the module
    print(combustion_obj)

    # Plot the results
    combustion_obj.plot_results()


def test_mass_simulator():
    """ test the mass simulator module with a simple test-case """

    # ----------- Generate objects:

    json_interpreter = generate_data_layer()
    mass_simulator_obj = System(json_interpreter.return_mass_simulator_table())

    # Print the total mass
    print("\nRockets Total Mass: {0} kgs".format(mass_simulator_obj.get_mass()))

    # Print the results
    print(mass_simulator_obj)


def test_trajectory():
    """ perform the test over the trajectory module """

    # -------------- Define parameters:

    # thrust & time
    delta_t = 0.01                                          # delta-time in seconds
    simulation_time = 45                                    # simulation-time
    n_points = int(simulation_time / delta_t) + 1           # total number of points
    thrust = []                                             # Initiate the thrust array
    time = []                                               # Initiate the time array
    burn_time = 3.5                                         # Burn time in seconds
    constant_thrust = 3000                                  # Thrust value in newtons
    for i in range(0, n_points):
        t = delta_t*i
        time.append(t)
        if t < burn_time:
            thrust.append(constant_thrust)
        else:
            thrust.append(0)

    # isp, area_ref, initial conditions
    isp = 245  # s
    initial_conditions = {'h0': 0, 'v0': 0, 'm0': 28}       # Initial mass in kg

    # ------------- Generate objects:

    json_interpreter = generate_data_layer()
    trajectory_data = json_interpreter.return_trajectory_table()
    density_obj = DensityLaw(trajectory_data['density'])

    drag_parameters = {'drag_coefficient': trajectory_data['drag']['Cd'],
                       'area_ref': trajectory_data['drag']['area_ref'],
                       'density': density_obj}

    drag_obj = DragCatalogue.return_drag_object(trajectory_data['drag']['type'], drag_parameters)
    trajectory_obj = TrajectoryObject(density_obj=density_obj, drag_obj=drag_obj)

    # -------------- Run simulation:

    trajectory_obj.run_simulation_on_trajectory(time=np.asarray(time),
                                                thrust=np.asarray(thrust),
                                                isp=isp,
                                                initial_conditions=initial_conditions)

    # Plot the results
    trajectory_obj.plot_results()

# ---------------------------- MAIN  ---------------------------------

if __name__ == '__main__':

    # Call on test_combustion method
    # test_combustion()
    # test_mass_simulator()
    # test_trajectory()
    test_combustion_onera_data()