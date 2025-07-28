# Core types and configuration options
add_subdirectory(core)

set(SHARE_SRC
  eamxx_config.cpp
  eamxx_session.cpp
  atm_process/atmosphere_process.cpp
  atm_process/atmosphere_process_hash.cpp
  atm_process/atmosphere_process_group.cpp
  atm_process/atmosphere_process_dag.cpp
  atm_process/atmosphere_diagnostic.cpp
  atm_process/IOPDataManager.cpp
  field/field_alloc_prop.cpp
  field/field_identifier.cpp
  field/field_header.cpp
  field/field_layout.cpp
  field/field_tracking.cpp
  field/field.cpp
  field/field_group.cpp
  field/field_manager.cpp
  field/field_sync.cpp
  grid/abstract_grid.cpp
  grid/grids_manager.cpp
  grid/grid_import_export.cpp
  grid/se_grid.cpp
  grid/point_grid.cpp
  grid/remap/abstract_remapper.cpp
  grid/remap/coarsening_remapper.cpp
  grid/remap/horiz_interp_remapper_base.cpp
  grid/remap/horiz_interp_remapper_data.cpp
  grid/remap/iop_remapper.cpp
  grid/remap/identity_remapper.cpp
  grid/remap/refining_remapper_p2p.cpp
  grid/remap/vertical_remapper.cpp
  property_checks/property_check.cpp
  property_checks/field_nan_check.cpp
  property_checks/field_within_interval_check.cpp
  property_checks/mass_and_energy_conservation_check.cpp
  util/eamxx_data_interpolation.cpp
  util/eamxx_fv_phys_rrtmgp_active_gases_workaround.cpp
  util/eamxx_time_interpolation.cpp
  util/eamxx_time_stamp.cpp
  util/eamxx_timing.cpp
  util/eamxx_utils.cpp
  util/eamxx_bfbhash.cpp
  util/eamxx_repro_sum_mod.F90
)

# Small/generic utilities
add_subdirectory(util)

# Basic physics constant/functions
add_subdirectory(physics)

# Field structures
add_subdirectory(field)

add_library(scream_share ${SHARE_SRC})
set_target_properties(scream_share PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/modules)
target_include_directories (scream_share PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/modules)

target_include_directories(scream_share PUBLIC
  ${SCREAM_SRC_DIR}
  ${SCREAM_BIN_DIR}/src
  ${CMAKE_CURRENT_SOURCE_DIR}
)

# Property checks structures
add_subdirectory(property_checks)

target_link_libraries(scream_share PUBLIC csm_share)

if (GPTL_PATH)
  target_include_directories(scream_share PUBLIC ${GPTL_PATH})
else()
  target_include_directories(scream_share PUBLIC ${INSTALL_SHAREDPATH}/include)
endif()

# Data managers structures
add_subdirectory(data_managers)

# Atm process
add_subdirectory(atm_process)

# Diagnostic classes
add_subdirectory(diagnostics)

# IO library
add_subdirectory(io)

# Algorithm library
add_subdirectory(algorithm)

# For bwd compatibility: a "catch" all interface lib of all the stuff in the share folder
add_library(scream_share INTERFACE)
target_link_libraries(scream_share INTERFACE
  eamxx_core
  eamxx_scorpio_interface
  eamxx_utils
  eamxx_field
  eamxx_grid
  eamxx_remap
  eamxx_property_checks
  eamxx_data_managers
  eamxx_atm_process
  eamxx_diagnostics
  eamxx_io
  eamxx_algorithm)

# The "build_cf_dictionary" target downloads an XML file containing valid field
# names and aliases based on the CF conventions, and transforms it into a YAML
# file that we use to validate our own fields against it. This target is only
# available if CMake finds wget on your system.
find_program(WGET wget)
if (NOT WGET STREQUAL "WGET-NOTFOUND")
  set(CF_CONVENTIONS_XML_FILE cf-standard-name-table.xml)
  set(CF_CONVENTIONS_URL "http://cfconventions.org/Data/cf-standard-names/76/src/${CF_CONVENTIONS_XML_FILE}")
  add_custom_target(build_cf_dictionary
    ${WGET} ${CF_CONVENTIONS_URL}
    COMMAND ${PROJECT_SOURCE_DIR}/scripts/cf-xml-to-yaml ${CF_CONVENTIONS_XML_FILE}
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/data
    COMMENT "Downloading CF conventions XML file and building YAML dictionary...")
endif()
