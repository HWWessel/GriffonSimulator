# testIntegration.py script is in charge of testing the integration of the different modules in
# order to allow the interaction between them.
# Author: Jose Felix Zapata Usandivaras
# Date: 2/1/2018.

# ---------------------- IMPORT MODULES -----------------------

from Initializer.Initializer import *                               # Import the Initializer object
from IntegrationModule.SimulationObject import SimulationObject     # Import the SimulationObject


# -------------------- FUNCTIONS DEFINITIONS ------------------


def generate_data_layer(data_file):
    """ generate_data_layer instantiates all of the objects that
    are required as data layer for the program to run properly.
    :param data_file: name of the json file to be used as reference for the simulation (str)
    :return JsonInterpreter """

    # Pass the file name
    data_directory = "C:/Users/Felix Zapata/Google Drive/" \
                     "MSc in AE ISAE - SUPAERO/Griffon Project/GriffonSimulator/data"

    # Return the output
    return JsonInterpreter(file_name="/".join([data_directory, data_file]))


def test_simulation_initializer():
    """ the method is in charge of generating the inputs for the initializer and
    instantiating it. """

    # ------------ Generate the data-layer:

    json_interpreter = generate_data_layer("Thermodynamic Data 36 bar OF 0,1 to 8,0 H2O2 87,5.json")

    # ---------- Pack the inputs:

    init_parameters = {
                        'combustion': {
                                       'geometric_params': {'type': OneCircularPort, 'L': 0.5, 'rintInitial': 0.04,
                                                            'rext0': 0.1},

                                       'nozzle_params': {'At': 0.000589, 'expansion': 5.7, 'lambda_e': 0.98,
                                                         'erosion': 0},

                                       'design_params': {'gamma': 1.27, 'p_chamber': 3600000, 'p_exit': 100000,
                                                         'c_star': 1500, 'isp': 230, 'thrust': 30000},
                                      },
                        'mass_simulator': {
                                            'tank_filling': 0.9
                                          }

                      }
    simulation_parameters = {
                              'combustion': {'ox_flow': 'TBD', 'injection_loss': 0.5, 'area_injection': 0.000105,
                                             'safety_thickness': 0.005, 'dt': 0.05},
                              'trajectory': {'initial_conditions': {'h0': 0, 'v0': 0, 'm0': 'TBD'},
                                             'simulation_time': 40}
                            }

    # -------------- Generate the initializer:

    init_obj = Initializer(init_parameters=init_parameters,
                           simulation_parameters=simulation_parameters,
                           json_interpreter=json_interpreter)

    # -------------- Generate the simulation object:

    simulation_object = SimulationObject(initializer_collection=init_obj)

    # --------------- Run the simulation:

    simulation_object.run_simulation_in_batch()

    # --------------- Plot the results

    # simulation_object.combustion_module.plot_results()
    # simulation_object.trajectory_module.plot_results()

    # print("Hello World"

    # Export results to csv files

    # data directory
    data_directory = "C:/Users/Felix Zapata/Google Drive/" \
                     "MSc in AE ISAE - SUPAERO/Griffon Project/GriffonSimulator/data/data_tests"

    file_name_expression = "Griffon Output Test {number}.csv"

    simulation_object.export_results_to_file(file_name_expression="/".join([data_directory,
                                                                            file_name_expression]))

    # Plot results

    simulation_object.results_collection.elements_list[0].combustion.plot_results()
    simulation_object.results_collection.elements_list[0].trajectory.plot_results()

# ----------------------------- MAIN ------------------------------

if __name__ == '__main__':

    # Call the test functions
    test_simulation_initializer()
