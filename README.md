# BlueSkyCMAQEmission
## Descriptions
The is a Python-based function to convert BlueSky pipeline outputs to 4D dimensional point-based fire emissions (TIME $\times$ X $\times$ Y $\times$ LAYERS). The point emission location is at the centroid of fire (provided by ```lat``` and ```lng``` attributes in BlueSky output).

BlueSky pipeline output should provide the intensity, timeprofile and vertical profile (plume bottom and plume top) to generate the emissions for CMAQ. The code inject the ***flaming*** phase of fire emissions into the vertical layers between plume bottom and plume top. The ***smodering*** phase and ***residual*** phase of fire emissions are injected into the surface layer.

Currently, the code only supports generating point-based fire emissions. If the BlueSky runs with the fire perimeter, the emission will treat the fire as a point source where the location is at the centroid of the fire perimeter. 

If you provide a METCRO3D which does not have time overlap with the BlueSky outputs, the code will generate a emission file with all zeros for all of the species.

## Usage
Use the function: ```write_fire_emissions(bsp_filename, metcros_filename, mechanism_species_mapping, output_filename)```
* ```bsp_filename```: the BlueSky pipeline outputs. BlueSky pipeline outputs provide the intensity, timeprofile, and vertical profile (plume bottom and plume top) of fire emissions.
* ```metcros_filename```: the METCRO3D file from MCIP. The METCRO3D file provides the vertical layer structure, the spatial domain range of CMAQ.
* ```mechanism_species_mapping```: the dictionary to map species from BlueSky pipeline output to species in CMAQ. The key in the dictionary is the species name in CMAQ (``str``). The value in the dictionary is a ``list`` of 2-element ``tuple``. The first element in the tuple is the coefficient (``float``) and the second element in the tuple is the species name in BlueSky (``str``). The mapping convention is shown below:
**Example**:
If there is a ``key-value`` as following:
``'CMAQ_species_1': [(coefficient_1, 'bsp_species_1'), (coefficient_2, 'bsp_species_2')]``
The mapping method for the ``key-value`` is:
``CMAQ_species_1 = (coefficient_1 *bsp_species_1) + (coefficient_2 * bsp_species_2)``
We provided two mapping mechanisms ``CB6_species_mapping`` and ``SAPRC07_species_mapping`` in ``Mechanism.py``.
* ```output_filename```: the name of the output NetCDF emission file.

## Code Structure
* EmissionGenerator.py: includes all functions for converting BleuSky outputs to CMAQ fire emissions.
* HeaderInfo.py: includes the variable list for emission files and descriptions of each variable. Users can use an anthropogenic emission to extract such information by using HeaderExtraction.py.
* HeaderExtraction.py: extracts the header information in the emission file (optional).
* Mechanism.py: definition of mapping BlueSky pollutants species to CMAQ under different chemical mechanisms.
* main.py: provides two examples for users to reference.
## Examples
We provide two examples for using the function. The examples are in ``main.py``.
* One is the prescribed fire simulation in Southeastern U.S. The fire lasts less than one day and one daily METCRO3D file can cover the period of prescribed fires.
* One is the wildfire fire simulation in Southeastern U.S. The fire lasts more than one day and need more than one daily METCRO3D file to cover the period of wildfire fires.

We also offer a CMAQ visualization that illustrates the appearance of the 4D fire emissions generated by the code:
<p align="center">
  <img src="https://github.com/zli867/BlueSkyCMAQEmission/blob/main/results/CMAQ_Briggs_smoke_0302.gif" />
</p>
