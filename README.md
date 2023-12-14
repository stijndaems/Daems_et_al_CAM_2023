# CAM metabolic proton homeostasis - scripts
This respository contains data and scripts used to study proton metabolism in mature CAM leaves
# Dependencies
libsbml version 5.12.1<br>
cobrapy version 0.5.4<br>
# Instructions
1. Install python version 2, libsbml and cobrapy<br>
2. Go to folder containing scripts<br>
3. Open a console (mac), terminal (linux) or command prompt (windows)<br>
4. Run scripts<br>
# Index
  * Model_construction: generate 12-phase model from PlantCoreMetabolism model, integrate 4 CAM phases of gas exchange, add metabolite accumulation constraints<br>
  * ME_refinement_and_simulations: refine the 12-phase model for a malic enzyme-type CAM leaf mesophyll cell, apply additional constraints as described in the paper, and run simulations<br>
  * Sensitivity_analysis: simulate the impact of PEPCK as main decarboxylase, varying Rubisco carboxylase/oxygenase (Vc/Vo) ratios, increased and decreased maintenance cost on key proton reactions + FVA<br>
  * Functions.py: a python module containing functions used in scripts Model_construction, ME_refinement_and_simulations and Sensitivity_analysis<br>
  * GasExchange.xlsx: experimental data of diel CO2 exchange reported in Phalaenopsis<br>
  * ProcessedData_v2.xlsx: experimental diel course metabolite data reported in Phalaenopsis used to apply metabolite accumulation constraints<br>
  * PlantCoreMetabolism_v2_0_0.xml: a proton balanced model of primary metabolism in plant cells<br>
