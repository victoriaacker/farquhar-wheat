# -*- coding: latin-1 -*-

import os
import pandas as pd

from farquharwheat import simulation, converter

from turgorgrowth import simulation as turgorgrowth_simulation, converter as turgorgrowth_converter, \
    postprocessing as turgorgrowth_postprocessing, tools as turgorgrowth_tools

"""
    main
    ~~~~

    An example to show how to:

        * initialize and run the model Farquhar-Wheat,
        * format the outputs of Farquhar-Wheat.

    You must first install :mod:`farquharwheat` and its dependencies
    before running this script with the command `python`.

    :copyright: Copyright 2014-2015 INRA-ECOSYS, see AUTHORS.
    :license: see LICENSE for details.

"""
# METEO_INPUTS_FILENAME = 'water_potential_forcings.csv'

INPUTS_ELEMENT_FILENAME = 'elements_inputs.csv'
INPUTS_AXIS_FILENAME = 'axes_inputs.csv'
OUTPUTS_FILENAME = 'outputs.csv'

OUTPUTS_PRECISION = 6

if __name__ == '__main__':

    # create a simulation and a converter
    simulation_ = simulation.Simulation()
    # read inputs from Pandas dataframe
    elements_inputs_df = pd.read_csv(INPUTS_ELEMENT_FILENAME, sep=",")
    axes_inputs_df = pd.read_csv(INPUTS_AXIS_FILENAME, sep=";")
    # convert the dataframe to simulation inputs format
    inputs = converter.from_dataframe(elements_inputs_df, axes_inputs_df)
    # initialize the simulation with the inputs
    simulation_.initialize(inputs)
    # run the simulation
    simulation_.run(Ta=18.8, RH=0.530000, Ur=2.200000, ambient_CO2=400, SRWC=80)
    # convert the outputs to Pandas dataframe
    outputs_df = converter.to_dataframe(simulation_.outputs)
    # write the dataframe to CSV
    outputs_df.to_csv(OUTPUTS_FILENAME, index=False, na_rep='NA', float_format='%.{}f'.format(OUTPUTS_PRECISION))
