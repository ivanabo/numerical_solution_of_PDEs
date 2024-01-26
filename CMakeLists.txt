add_executable("MF2" MF2.cc operator_MF2.hh bctype_MF2.hh combo_driver.hh time_operator.hh space_operator.hh permeability_generatorMF2.hh slvel.hh)
target_link_dune_default_libraries("MF2")
