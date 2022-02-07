# Copy h5 files to binary dir
file(GLOB_RECURSE all_h5_ref_files RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.ref.h5)
foreach(file ${all_h5_ref_files})
  configure_file(${file} ${file} COPYONLY)
endforeach()

# List of all tests
set(all_tests
  basic/ct_seg_single_bath
  basic/ct_seg_single_bath_move_move
  basic/ct_seg_single_bath_improved_est
  basic/ct_seg_single_bath_retint
  basic/ct_seg_single_bath_spinspin
  basic/ct_seg_single_bath_spinspin_move_move
  basic/ct_seg_multiorb
  offdiag/offdiag_ctseg
  offdiag/offdiag_ctseg_4_colors
  spinboson/spinspin_subohmic
  vertex/ct_seg_single_bath_vertex
)

foreach(test ${all_tests})
  get_filename_component(test_name ${test} NAME_WE)
  get_filename_component(test_dir ${test} DIRECTORY)
  add_test(NAME Py_${test_name} COMMAND ${TRIQS_PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/${test_dir}/${test_name}.py WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${test_dir})
  set_property(TEST Py_${test_name} APPEND PROPERTY ENVIRONMENT PYTHONPATH=${PROJECT_BINARY_DIR}/python:$ENV{PYTHONPATH} ${SANITIZER_RT_PRELOAD})
endforeach()
