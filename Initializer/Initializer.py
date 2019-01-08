# Initializer.py contains the Initializer class implementation in order to achieve the initialization of the different
# modules of the Griffon program, that is Combustion and MassSimulator module in principle.
# Author: Jose Felix Zapata Usandivaras
# Date: 12/12/2018


# ------------------- IMPORT MODULES ----------------------------

from math import pi, sqrt                                   # Import math functions
from CombustionModule.Geometries import *                   # Import the Geometries
from CombustionModule.Nozzle import *                       # Import the nozzles
from DataLayer.JsonInterpreter import JsonInterpreter       # Import the JsonInterpreter class
from TrajectoryModule.Drag import *                         # Import the Drag Module (DensityLaw in it)
from Libraries.Collections import Collections               # Import Collections class
# ------------------ FUNCTIONS/STATIC CLASSES -------------------


class InitializerCalculator:
    """
    InitializerCalculator is an auxiliary class which helps decoupling the different
    calculations present in the code and handle then separately to ease debugging.
    It acts as a library of methods/function. All the presented methods are static and it
    does not have and initializer.
    """

    @staticmethod
    def calculate_fuel_mass(geometry_obj, fuel_density):
        """ calculate_fuel_mass determines the mass of fuel contained in the combustion chamber.
        :param geometry_obj: Geometry instance used to calculate initial mass
        :param fuel_density: density of the fuel in the combustion chamber [kg/m^3] """

        # Check the inputs
        assert isinstance(geometry_obj, Geometry), "Please insert a valid Geometry type. \n"

        # Return the output
        return geometry_obj.get_fuel_mass(fuel_density)

    @staticmethod
    def calculate_oxidiser_flow(oxidant_pressure, chamber_pressure, injection_loss, area_injection, ox_density):
        """ calculate_oxidiser_flow determines the flow of oxidizer based on chamber and oxidant_tank pressure.
        :param oxidant_pressure: pressure in [Pa] at the oxidant tank.
        :param chamber_pressure: pressure in [Pa] at the combustion chamber.
        :param injection_loss: coefficient for injection 0 < K* < 1, as a proportion of Delta P
        :param area_injection: area of injection [m^2]
        :param ox_density: density of the oxidiser [kg/m^3] """

        # Return the oxidizer flow
        return area_injection * sqrt(2 * ox_density * (oxidant_pressure - chamber_pressure) * (1 - injection_loss))

    @staticmethod
    def calculate_oxidiser_mass(tank_radius, tank_height, oxidiser_density, tank_filling):
        """
        calculate_oxidiser_mass determines the mass of oxidizer for a complete full tank
        :param tank_radius: radius of the tank
        :param tank_height: height of the tank
        :param oxidiser_density: density of the oxidizer
        :param tank_filling: % of the tank that is filled
        :return: total o
        """

        # TODO: this issue should be contemplated in the part.py module with an appropiate class.
        return pi * tank_radius * tank_radius * tank_height * tank_filling * oxidiser_density

    @staticmethod
    def calculate_oxidiser_density(ox_density_relative, ox_purity):
        """
        calculate_oxidiser_density determines the oxidiser density based on purity and actual density
        of the pure component.
        :param ox_density_relative: relative density of the pure oxidizer (relative to water)
        :param ox_purity: purity of the oxidizer (mixed with water).
        :return: density of the oxidizer [kg/m^3]
        """
        water_density = 1000            # Define the water density
        return water_density / (1 - ox_purity * (1 - 1 / ox_density_relative))

# --------------------- CLASS DEFINITIONS -----------------------


class Initializer:
    """
    The Initializer class pretends to generate the inputs required by the combustion
    and the MassSimulator module.

    Attributes:
        1. combustion_parameters: store the calculated combustion_parameters dictionary
        2. mass_simulator_parameters: store the calculated mass_simulator_parameters
        3. trajectory_parameters: store the calculated trajectory_parameters
        4. simulation_parameters: store the simulation parameters in order to command the simulation run
    """

    def __init__(self, init_parameters, simulation_parameters, json_interpreter):
        """
        class initializer
        :param init_parameters: list of parameters defined to initialize the complete
        object as a hole.
        :param simulation_parameters: parameters associated to the run of the simulation
        :param json_interpreter: JsonInterpreter instance (singleton)
        """

        # Check the inputs
        assert isinstance(init_parameters, dict), "The parameters of the Simulation must be passed as a dict \"" \
                                                        "instance \n"
        assert json_interpreter == JsonInterpreter.instance, "Please insert a valid JsonInterpreter instance. \n"

        # Store the input as an attribute
        self.combustion_parameters = self._generate_combustion_parameters(init_parameters['combustion'],
                                                                          json_interpreter)
        self.mass_simulator_parameters = self._generate_mass_simulator_parameters(init_parameters['mass_simulator'],
                                                                                  json_interpreter)
        self.trajectory_parameters = self._generate_trajectory_parameters(json_interpreter)

        self.simulation_parameters = self._state_simulation_parameters(simulation_parameters, json_interpreter)

    def __str__(self):
        """ return a string representation of the object based on the geometry, nozzle and simulation parameters
        implemented. """

        # Get the geometry and nozzle strings
        geometry_obj = self.combustion_parameters['geometry_object']
        nozzle_obj = self.combustion_parameters['nozzle_object']

        # Introduction string
        intro_str = "\nInitializing simulation w/parameters:"

        # Generate a string for the geometry and nozzle
        inputs_str = "\n{geom} \n {nozzle}".format(geom=geometry_obj, nozzle=nozzle_obj)

        # Simulation parameters string
        data_str = "\nSimulation Parameters: \n\t" + \
                   "\t\n\t".join(("{name}: {value}".format(name=variable, value=val) for variable, val
                                  in self.simulation_parameters.items()))

        return "\n".join([intro_str, inputs_str, data_str])

    @staticmethod
    def _generate_combustion_parameters(combustion_init_dict, json_interpreter):
        """
        _generate_combustion_parameters performs all of the necessary calculations (if there are )required to set the
        inputs for the combustion module
        :param combustion_init_dict: dictionary containing the data required to prepare the inputs for the
        combustion module.
        :param json_interpreter: JsonInterpreter instance
        :return: combustion_parameters dictionary
        """

        # ------------- Generate the objects:

        # Geometry
        geometry_obj = combustion_init_dict['geometric_params'].pop('type')
        geometry_obj = geometry_obj(**combustion_init_dict['geometric_params'])

        # Nozzle
        nozzle_obj = Nozzle(**combustion_init_dict['nozzle_params'])
        nozzle_obj.set_design(**combustion_init_dict['design_params'])

        # Return output
        return {'json_interpreter': json_interpreter, 'geometry_object': geometry_obj, 'nozzle_object': nozzle_obj}

    def _generate_mass_simulator_parameters(self, mass_simulator_init_dict, json_interpreter):
        """
        _generate_mass_simulator_parameters performs all of the necessary calculations (if there are) required to set
        the inputs for the mass simulation module
        :param mass_simulator_init_dict: initializer dictionary for the mass_simulator
        :param json_interpreter: JsonInterpreter instance
        :return: mass estimation parameters dictionary
        """

        # Extract the system dictionary with JsonInterpreter
        system_dict = json_interpreter.return_mass_simulator_table()

        # Get combustion table
        combustion_table = json_interpreter.return_combustion_table()

        # Calculate propellant & oxidiser mass
        chamber_length = self.combustion_parameters['geometry_object'].length
        chamber_radius = self.combustion_parameters['geometry_object'].r_ext
        propellant_mass = InitializerCalculator.calculate_fuel_mass(self.combustion_parameters['geometry_object'],
                                                                    combustion_table['rho_fuel'])

        ox_density = InitializerCalculator.calculate_oxidiser_density(combustion_table['rho_ox_pure'],
                                                                      combustion_table['ox_purity'])

        tank_r = system_dict['subsystems']['oxidiser']['parts']['oxidant_tank']['radius']
        tank_h = system_dict['subsystems']['oxidiser']['parts']['oxidant_tank']['height']
        tank_f = mass_simulator_init_dict['tank_filling']
        ox_mass = InitializerCalculator.calculate_oxidiser_mass(tank_radius=tank_r,
                                                                tank_height=tank_h,
                                                                tank_filling=tank_f,
                                                                oxidiser_density=ox_density)

        # Set the system_dict to its proper values
        system_dict['subsystems']['combustion']['parts']['chamber']['height'] = chamber_length
        system_dict['subsystems']['combustion']['parts']['chamber']['radius'] = chamber_radius
        system_dict['subsystems']['combustion']['parts']['chamber']['propellant_mass'] = propellant_mass
        system_dict['subsystems']['oxidiser']['parts']['oxidant_tank']['propellant_mass'] = ox_mass

        # Return the output
        return {'system_dict': system_dict}

    @staticmethod
    def _generate_trajectory_parameters(json_interpreter):
        """
        _generate_trajectory_parameters generates the inputs for the TrajectoryModule
        :param json_interpreter: JsonInterpreter instance (singleton)
        :return: TrajectoryModule inputs dictionary
        """

        # TODO: review if possible to implement calculation of aerodynamic area of rocket based on system dict

        # Extract the trajectory data and generate the objects
        trajectory_data = json_interpreter.return_trajectory_table()
        density_obj = DensityLaw(trajectory_data['density'])

        drag_parameters = {'drag_coefficient': trajectory_data['drag']['Cd'],
                           'area_ref': trajectory_data['drag']['area_ref'],
                           'density': density_obj}
        drag_obj = DragCatalogue.return_drag_object(trajectory_data['drag']['type'], drag_parameters)

        # Return the dictionary
        return {'density_obj': density_obj, 'drag_obj': drag_obj}

    def _state_simulation_parameters(self, simulation_parameters, json_interpreter):
        """
        _state_simulation_parameters helps determine the simulation parameters that are linked and calculate
        them in order to run the simulation. In the present case we are talking about oxidizer flow.
        :param simulation_parameters: simulation_parameters dictionary
        :param json_interpreter: JsonInterpreter instance (singleton)
        :return: modified simulation_parameters dictionary
        """

        # ------------ Recalculate oxidizer flow

        # Extract the data
        system_dict = self.mass_simulator_parameters['system_dict']
        combustion_table = json_interpreter.return_combustion_table()

        # Oxidant pressure
        ox_pressure = system_dict['subsystems']['oxidiser']['parts']['oxidant_tank']['pressure']

        # Chamber pressure
        chamber_pressure_mass = system_dict['subsystems']['combustion']['parts']['chamber']['pressure']
        chamber_pressure_combustion = combustion_table['P_chamber_bar'] * (10 ** 5)

        assert chamber_pressure_combustion == chamber_pressure_mass, "Chamber pressure in both systems are not the " \
                                                                     "same. \n"

        # TODO: review if there is a better way to store injection loss and area of injection (combustion table/other)
        # Injection loss
        injection_loss = simulation_parameters['combustion'].pop('injection_loss')

        # Area injection
        area_injection = simulation_parameters['combustion'].pop('area_injection')

        # Oxidiser Density
        ox_density = InitializerCalculator.calculate_oxidiser_density(combustion_table['rho_ox_pure'],
                                                                      combustion_table['ox_purity'])

        # Calculate ox_flow
        ox_flow = InitializerCalculator.calculate_oxidiser_flow(oxidant_pressure=ox_pressure,
                                                                chamber_pressure=chamber_pressure_combustion,
                                                                injection_loss=injection_loss,
                                                                area_injection=area_injection,
                                                                ox_density=ox_density)

        # Allocate the value
        simulation_parameters['combustion']['ox_flow'] = ox_flow

        # Return the dictionary
        return simulation_parameters


class InitializerCollection(Collections):
    """
    InitializerCollection class works as a container of Initializer
    required to perform general operations to multiple Initializer
    Attributes:
        0. elements: list containing the constellation
    """

    def __init__(self, initializer_obj):
        """
        class initializer
        :param initializer_obj: Initializer instance, iterable of Initializer, or empty list
        """
        # Call superclass method
        super().__init__(initializer_obj, Initializer)

    def add_element(self, element, *args):
        """ Override method from parent class"""
        super().add_element(element, Initializer)
