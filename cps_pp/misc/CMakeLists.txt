set( RPCGENVML_SRC 
  vml_encode/rpcgenvml/rpc_clntout.c
  vml_encode/rpcgenvml/rpc_cmsg.c
  vml_encode/rpcgenvml/rpc_common.c
  vml_encode/rpcgenvml/rpc_cout.c
  vml_encode/rpcgenvml/rpc_dtable.c
  vml_encode/rpcgenvml/rpc_hout.c
  vml_encode/rpcgenvml/rpc_main.c
  vml_encode/rpcgenvml/rpc_parse.c
  vml_encode/rpcgenvml/rpc_pktout.c
  vml_encode/rpcgenvml/rpc_prot.c
  vml_encode/rpcgenvml/rpc_sample.c
  vml_encode/rpcgenvml/rpc_scan.c
  vml_encode/rpcgenvml/rpc_svcout.c
  vml_encode/rpcgenvml/rpc_tblout.c
  vml_encode/rpcgenvml/rpc_util.c
  vml_encode/rpcgenvml/stpcpy.c )

SET(CMAKE_C_COMPILER gcc)
#set(CMAKE_EXE_LINKER_FLAGS -openmp)

add_executable( rpcgenvml ${RPCGENVML_SRC} )

set( ARG_VML_SRC
  ape_smear_arg
  array_arg
  bfm_arg
  cg_arg
  dens_arg
  do_arg
  eig_arg
  eigcg_arg
  enum
  eq_state_arg
  fix_gauge_arg
  ghb_arg
  hmc_arg
  hmd_arg
  hyp_smear_arg
  inst_arg
  int_arg
  kinetic_smear_arg
  lanc_arg
  mdwf_arg
  meas_arg
  mom_arg
  noise_arg
  nuc3pt_arg
  pbp_arg
  pot_arg
  qpropw_arg
  remez_arg
  s_spect_arg
  w_spect_arg )


FOREACH(item ${ARG_VML_SRC})
  add_custom_command(
    OUTPUT "${CMAKE_BINARY_DIR}/generated/include/${item}.h"
    OUTPUT "${CMAKE_BINARY_DIR}/generated/src/${item}.C"
    COMMAND cd ${CMAKE_SOURCE_DIR}/include/input && ${CMAKE_CURRENT_BINARY_DIR}/rpcgenvml -h -C -x -u util/vml/ -v alg/  ${item}.x >${CMAKE_BINARY_DIR}/generated/include/alg/${item}.h
    COMMAND cd ${CMAKE_SOURCE_DIR}/include/input && ${CMAKE_CURRENT_BINARY_DIR}/rpcgenvml -c -C -x -u util/vml/ -v alg/  ${item}.x >${CMAKE_BINARY_DIR}/generated/src/${item}.C
    DEPENDS "${CMAKE_SOURCE_DIR}/include/input/${item}.x" rpcgenvml
  )
  LIST( APPEND X_LIST "${CMAKE_BINARY_DIR}/generated/include/${item}.h" )
  LIST( APPEND X_LIST "${CMAKE_BINARY_DIR}/generated/src/${item}.C" )
ENDFOREACH(item ${ARG_VML_SRC})

add_custom_target(
  generate-foo
  DEPENDS ${X_LIST}
)


