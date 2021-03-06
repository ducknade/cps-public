set(CPS_OBJS
  alg/alg_actiondensity/alg_actiondensity.C
  alg/alg_arg/ape_smear_arg_init.C
  alg/alg_arg/array_arg_init.C
  alg/alg_arg/cg_arg_init.C
  alg/alg_arg/do_arg_init.C
  alg/alg_arg/eig_arg_init.C
  alg/alg_arg/eigcg_arg_init.C
  alg/alg_arg/enum_int.C
  alg/alg_arg/hyp_smear_arg_init.C
  alg/alg_arg/int_arg_resize.C
  alg/alg_arg/kinetic_smear_arg_init.C
  alg/alg_arg/nuc3pt_arg_init.C
  alg/alg_arg/qpropw_arg_init.C
  alg/alg_arg/remez_arg_init.C
  alg/alg_arg/threept_arg.C
  alg/alg_arg/threept_arg_init.C
  alg/alg_arg/w_spect_arg_init.C
  alg/alg_base/alg_base.C
  alg/alg_dens/alg_dens.C
  alg/alg_eig/alg_eig.C
  alg/alg_eq_state/alg_eq_state.C
  alg/alg_fix_gauge/alg_fix_gauge.C
  alg/alg_fourier_prop/FourMom.C
  alg/alg_fourier_prop/alg_fourier_prop.C
  alg/alg_hmd/alg_hmc_phi.C
  alg/alg_hmd/alg_hmc_qpq.C
  alg/alg_hmd/alg_hmc_rhmc.C
  alg/alg_hmd/alg_hmd.C
  alg/alg_hmd/alg_hmd_r.C
  alg/alg_hmd/alg_hmd_r2.C
  alg/alg_hq_pot/alg_hq_pot.C
  alg/alg_inst/alg_inst.C
  alg/alg_int/alg_action.C
  alg/alg_int/alg_action_bilinear.C
  alg/alg_int/alg_action_boson.C
  alg/alg_int/alg_action_fermion.C
  alg/alg_int/alg_action_gauge.C
  alg/alg_int/alg_action_quotient.C
  alg/alg_int/alg_action_rational.C
  alg/alg_int/alg_action_rational_quotient.C
  alg/alg_int/alg_action_rational_split.C
  alg/alg_int/alg_hamiltonian.C
  alg/alg_int/alg_hmc.C
  alg/alg_int/alg_int.C
  alg/alg_int/alg_int_ab.C
  alg/alg_int/alg_int_campostrini.C
  alg/alg_int/alg_int_force_grad.C
  alg/alg_int/alg_int_leap.C
  alg/alg_int/alg_int_omelyan.C
  alg/alg_int/alg_int_omelyan_44.C
  alg/alg_int/alg_int_omelyan_45.C
  alg/alg_int/alg_int_sum.C
  alg/alg_int/alg_momentum.C
  alg/alg_int/force_arg.C
  alg/alg_meas/alg_meas.C
  alg/alg_mom/alg_mom.C
  alg/alg_noise/alg_noise.C
  alg/alg_nuc3pt/CorrFunc.C
  alg/alg_nuc3pt/DTerms.C
  alg/alg_nuc3pt/Derivative.C
  alg/alg_nuc3pt/Gamma.C
  alg/alg_nuc3pt/Meson.C
  alg/alg_nuc3pt/Nuc2pt.C
  alg/alg_nuc3pt/Nuc3pt.C
  alg/alg_nuc3pt/Nuc3ptGammaR.C
  alg/alg_nuc3pt/alg_nuc3pt.C
  alg/alg_nuc3pt/qio_prop.C
  alg/alg_overrelax/comsrc/alg_overrelax.C
  alg/alg_overrelax/noarch/alg_overrelax_kern.C
  alg/alg_pbp/alg_pbp.C
  alg/alg_plaq/alg_plaq.C
  alg/alg_pot/alg_pot.C
  alg/alg_pot2/alg_pot2.C
  alg/alg_qpropw/Diquark.C
  alg/alg_qpropw/FermionVector.C
  alg/alg_qpropw/QPropW.C
  alg/alg_qpropw/QPropWFactory.C
  alg/alg_qpropw/SpinMatrix.C
  alg/alg_qpropw/WilsonMatrix.C
  alg/alg_qpropw/eig_QPropW.C
  alg/alg_qpropw/eq_mult.C
  alg/alg_rect/alg_rect.C
  alg/alg_remez/alg_remez.C
  alg/alg_rnd_gauge/alg_rnd_gauge.C
  alg/alg_s_spect/alg_s_spect.C
  alg/alg_s_spect/aots_s.C
  alg/alg_s_spect/hadron_prop_s.C
  alg/alg_s_spect/meson_prop_s.C
  alg/alg_s_spect/mom_meson_p_s.C
  alg/alg_s_spect/nlocal_prop_s.C
  alg/alg_s_spect/nlocal_propmes_s.C
  alg/alg_s_spect/nucl_prop_s.C
  alg/alg_s_spect/quark_prop_s.C
  alg/alg_smear/alg_smear.C
  alg/alg_smear/alg_smear_oleg.C
  alg/alg_smear2/alg_smear2.C
  alg/alg_tcharge/alg_tcharge.C
  alg/alg_w_spect/alg_w_spc_ex.C
  alg/alg_w_spect/alg_w_spect.C
  alg/alg_w_spect/w_axialcurr.C
  alg/alg_w_spect/w_baryon.C
  alg/alg_w_spect/w_ext_mesonBE.C
  alg/alg_w_spect/w_ext_mesons.C
  alg/alg_w_spect/w_ferm_vec.C
  alg/alg_w_spect/w_field.C
  alg/alg_w_spect/w_fuzzing.C
  alg/alg_w_spect/w_gamma_mat.C
  alg/alg_w_spect/w_ginfo.C
  alg/alg_w_spect/w_hyper_rect.C
  alg/alg_w_spect/w_mesons.C
  alg/alg_w_spect/w_momenta.C
  alg/alg_w_spect/w_quark.C
  alg/alg_w_spect/w_quark_d.C
  alg/alg_wilsonflow/alg_wilsonflow.C
  alg/alg_wilsonflow/su3generator.C
  alg/alg_wline/alg_wline.C
  alg/ghb/comsrc/ghb.C
  alg/ghb/noarch/ghb_krn.C
  alg/ghb/noarch/ghb_sup.C
  alg/plaq/plaq.C
  comms/noarch/glb_sim/check_sum_sim.C
  comms/noarch/glb_sim/glb_min_max_sim.C
  comms/noarch/glb_sim/glb_sum_dir_sim.C
  comms/noarch/glb_sim/glb_sum_five_sim.C
  comms/noarch/glb_sim/glb_sum_matrix_dir_sim.C
  comms/noarch/glb_sim/glb_sum_multi_dir_sim.C
  comms/noarch/glb_sim/glb_sum_sim.C
  comms/noarch/glb_sim/slice_sum_sim.C
  comms/noarch/scu_sim/get_data_sim.C
  comms/qmp/cbuf/cbuf.C
  comms/qmp/cbuf/cram.C
  comms/qmp/glb_cpp/glb_min_max.C
  comms/qmp/glb_cpp/glb_sum.C
  comms/qmp/glb_cpp/glb_sum_dir.C
  comms/qmp/glb_cpp/glb_sum_five.C
  comms/qmp/glb_cpp/glb_sum_internal.C
  comms/qmp/glb_cpp/glb_sum_matrix_dir.C
  comms/qmp/glb_cpp/glb_sum_multi_dir.C
  comms/qmp/glb_cpp/slice_sum.C
  comms/qmp/scu/check_sum.C
  comms/qmp/scu/get_data.C
  comms/qmp/scu/glb_sum2.C
  comms/qmp/sysfunc.C
  mem/p2v/p2v.C
  util/alloc/comsrc/amalloc.C
  util/alloc/comsrc/smalloc_common.C
  util/alloc/noarch/fmalloc.C
  util/alloc/noarch/pmalloc.C
  util/alloc/noarch/smalloc.C
  util/checksum/comsrc/checksum.C
  util/command_line/command_line.C
  util/data_shift/comsrc/data_shift.C
  util/data_shift/noarch/shift_util.C
  util/dense_matrix/dense_matrix.C
  util/dirac_op/d_op_asqtad/noarch/asqtad_dirac.C
  util/dirac_op/d_op_asqtad/noarch/d_op_asqtad.C
  util/dirac_op/d_op_asqtad/noarch/dirac_serial_cpp.C
  util/dirac_op/d_op_base/comsrc/dirac_op_base.C
  util/dirac_op/d_op_base/comsrc/inv_cg.C
  util/dirac_op/d_op_base/comsrc/inv_eig_cg.C
  util/dirac_op/d_op_base/comsrc/jacobi.C
  util/dirac_op/d_op_base/comsrc/minresext.C
  util/dirac_op/d_op_base/comsrc/ritz.C
  util/dirac_op/d_op_base/noarch/inv_cg.C
  util/dirac_op/d_op_base/noarch/inv_cg_blas.C
  util/dirac_op/d_op_base/noarch/minvcg.C
  util/dirac_op/d_op_clover/noarch/clover.C
  util/dirac_op/d_op_clover/noarch/d_op_clover.C
  util/dirac_op/d_op_clover/noarch/d_op_clover_supp.C
  util/dirac_op/d_op_dwf/comsrc/d_op_dwf_cg_dwf.C
  util/dirac_op/d_op_dwf/d_op_dwf.C
  util/dirac_op/d_op_dwf/eig_MatInv.C
  util/dirac_op/d_op_dwf/noarch/dwf_dslash.C
  util/dirac_op/d_op_dwf/noarch/dwf_dslash_4.C
  util/dirac_op/d_op_dwf/noarch/dwf_dslash_5_plus.C
  util/dirac_op/d_op_dwf/noarch/dwf_dslash_5_plus_omp.C
  util/dirac_op/d_op_dwf/noarch/dwf_end.C
  util/dirac_op/d_op_dwf/noarch/dwf_init.C
  util/dirac_op/d_op_dwf/noarch/dwf_m.C
  util/dirac_op/d_op_dwf/noarch/dwf_mdag.C
  util/dirac_op/d_op_dwf/noarch/dwf_mdagm.C
  util/dirac_op/d_op_mdwf/d_op_mdwf.C
  util/dirac_op/d_op_p4/noarch/d_op_p4.C
  util/dirac_op/d_op_p4/noarch/p4_dirac.C
  util/dirac_op/d_op_stag/noarch/d_op_stag.C
  util/dirac_op/d_op_stag/noarch/dirac.C
  util/dirac_op/d_op_stag_types/d_op_stag_types.C
  util/dirac_op/d_op_stag_types/eigen_stag.C
  util/dirac_op/d_op_wilson/noarch/d_op_wilson.C
  util/dirac_op/d_op_wilson/noarch/wilson_dslash.C
  util/dirac_op/d_op_wilson/noarch/wilson_end.C
  util/dirac_op/d_op_wilson/noarch/wilson_init.C
  util/dirac_op/d_op_wilson/noarch/wilson_m.C
  util/dirac_op/d_op_wilson/noarch/wilson_mdag.C
  util/dirac_op/d_op_wilson/noarch/wilson_mdagm.C
  util/dirac_op/d_op_wilsonTm/d_op_wilsonTm.C
  util/dirac_op/d_op_wilsonTm/d_op_wilsonTm_bfm.C
  util/dirac_op/d_op_wilson_types/bicgstab.C
  util/dirac_op/d_op_wilson_types/d_op_wilson_types.C
  util/dirac_op/d_op_wilson_types/eigen_wilson.C
  util/error/error.C
  util/error/noarch/error_hdwcheck.C
  util/gauge_field/gauge_field.C
  util/gauge_field/gauge_field_functions.C
  util/gjp/gjp.C
  util/gjp/qmp/start_end.C
  util/lattice/convert/common.C
  util/lattice/convert/convert.C
  util/lattice/convert/convert_func.C
  util/lattice/f_asqtad/comsrc/Fforce.C
  util/lattice/f_asqtad/comsrc/Fforce_utils.C
  util/lattice/f_asqtad/comsrc/RHMC_Fforce.C
  util/lattice/f_asqtad/comsrc/f_asqtad.C
  util/lattice/f_asqtad/noarch/force_product_sum.C
  util/lattice/f_bfm/bfm_eigcg_linalg.C
  util/lattice/f_bfm/f_bfm.C
  util/lattice/f_clover/f_clover.C
  util/lattice/f_dwf/f_dwf.C
  util/lattice/f_dwf_base/f_dwf_base.C
  util/lattice/f_dwf_base/noarch/f_dwf_base_force.C
  util/lattice/f_mdwf/f_mdwf.C
  util/lattice/f_none/f_none.C
  util/lattice/f_p4/comsrc/Fforce_p4.C
  util/lattice/f_p4/comsrc/Fforce_utils_p4.C
  util/lattice/f_p4/comsrc/RHMC_Fforce_p4.C
  util/lattice/f_p4/comsrc/f_p4.C
  util/lattice/f_p4/noarch/force_product_sum_p4.C
  util/lattice/f_p4/noarch/p4_Fforce.C
  util/lattice/f_smear/f_smear.C
  util/lattice/f_stag/f_stag.C
  util/lattice/f_stag_types/f_stag_t.C
  util/lattice/f_wilson/f_wilson.C
  util/lattice/f_wilsonTm/f_wilsonTm.C
  util/lattice/f_wilson_types/comsrc/f_wilson_t.C
  util/lattice/f_wilson_types/comsrc/sigmaproj_tr.C
  util/lattice/f_wilson_types/fforce_wilson_type.C
  util/lattice/f_wilson_types/noarch/gamma_5.C
  util/lattice/f_wilson_types/noarch/sproj_tr.C
  util/lattice/fix_gauge/fix_gauge.C
  util/lattice/g_impr_OLSym/g_impr_OLSym.C
  util/lattice/g_impr_rect/g_impr_rect.C
  util/lattice/g_impr_rect/noarch/g_impr_rect_force.C
  util/lattice/g_none/g_none.C
  util/lattice/g_power_plaq/g_power_plaq.C
  util/lattice/g_power_rect/g_power_rect.C
  util/lattice/g_tadpole_rect/g_tadpole_rect.C
  util/lattice/g_wilson/g_wilson.C
  util/lattice/g_wilson/noarch/g_wilson_force.C
  util/lattice/lattice_base/lat_cont.C
  util/lattice/lattice_base/lattice_base.C
  util/lattice/lattice_base/link_buffer.C
  util/lattice/lattice_types/lattice_na.C
  util/lattice/lattice_types/lattice_nb.C
  util/lattice/lattice_types/lattice_nc.C
  util/lattice/lattice_types/lattice_nd.C
  util/lattice/lattice_types/lattice_nm.C
  util/lattice/lattice_types/lattice_nn.C
  util/lattice/lattice_types/lattice_np.C
  util/lattice/lattice_types/lattice_ns.C
  util/lattice/lattice_types/lattice_nw.C
  util/lattice/lattice_types/lattice_oa.C
  util/lattice/lattice_types/lattice_op.C
  util/lattice/lattice_types/lattice_pa.C
  util/lattice/lattice_types/lattice_pc.C
  util/lattice/lattice_types/lattice_pd.C
  util/lattice/lattice_types/lattice_pm.C
  util/lattice/lattice_types/lattice_pn.C
  util/lattice/lattice_types/lattice_ps.C
  util/lattice/lattice_types/lattice_pw.C
  util/lattice/lattice_types/lattice_ra.C
  util/lattice/lattice_types/lattice_rb.C
  util/lattice/lattice_types/lattice_rc.C
  util/lattice/lattice_types/lattice_rd.C
  util/lattice/lattice_types/lattice_rm.C
  util/lattice/lattice_types/lattice_rn.C
  util/lattice/lattice_types/lattice_rp.C
  util/lattice/lattice_types/lattice_rs.C
  util/lattice/lattice_types/lattice_rw.C
  util/lattice/lattice_types/lattice_sc.C
  util/lattice/lattice_types/lattice_sd.C
  util/lattice/lattice_types/lattice_sn.C
  util/lattice/lattice_types/lattice_ss.C
  util/lattice/lattice_types/lattice_sw.C
  util/lattice/lattice_types/lattice_ta.C
  util/lattice/lattice_types/lattice_tn.C
  util/lattice/lattice_types/lattice_tp.C
  util/lattice/lattice_types/lattice_va.C
  util/lattice/lattice_types/lattice_vc.C
  util/lattice/lattice_types/lattice_vd.C
  util/lattice/lattice_types/lattice_vn.C
  util/lattice/lattice_types/lattice_vs.C
  util/lattice/lattice_types/lattice_vw.C
  util/lattice/lattice_types/lattice_wa.C
  util/lattice/lattice_types/lattice_wc.C
  util/lattice/lattice_types/lattice_wd.C
  util/lattice/lattice_types/lattice_wm.C
  util/lattice/lattice_types/lattice_wn.C
  util/lattice/lattice_types/lattice_wp.C
  util/lattice/lattice_types/lattice_ws.C
  util/lattice/lattice_types/lattice_ww.C
  util/libmarshall/vml.C
  util/libmarshall/vml_encoder.C
  util/libmarshall/xdr_encoder.C
  util/libmarshall/xml_encoder.C
  util/mom/mom.C
  util/momentum/momentum.C
  util/parallel_transport/pt_asqtad/noarch/pt_asqtad.C
  util/parallel_transport/pt_base/bagel_c/copy_gauge.C
  util/parallel_transport/pt_base/bagel_c/copy_gauge_s.C
  util/parallel_transport/pt_base/bagel_c/copy_matrix.C
  util/parallel_transport/pt_base/bagel_c/copy_matrix_s.C
  util/parallel_transport/pt_base/bagel_c/cross_lin.C
  util/parallel_transport/pt_base/bagel_c/cross_lin_s.C
  util/parallel_transport/pt_base/bagel_c/cross_look.C
  util/parallel_transport/pt_base/bagel_c/cross_look_s.C
  util/parallel_transport/pt_base/bagel_c/cross_over_lin.C
  util/parallel_transport/pt_base/bagel_c/cross_over_lin_s.C
  util/parallel_transport/pt_base/bagel_c/cross_over_look.C
  util/parallel_transport/pt_base/bagel_c/cross_over_look_s.C
  util/parallel_transport/pt_base/noarch/pt.C
  util/parallel_transport/pt_base/pt_base.C
  util/parallel_transport/pt_base/qmp/pt_alloc.C
  util/parallel_transport/pt_base/qmp/pt_amalloc.C
  util/parallel_transport/pt_base/qmp/pt_c_binding.C
  util/parallel_transport/pt_base/qmp/pt_cpp.C
  util/parallel_transport/pt_base/qmp/pt_force.C
  util/parallel_transport/pt_base/qmp/pt_init.C
  util/parallel_transport/pt_base/qmp/pt_init_g.C
  util/parallel_transport/pt_base/qmp/pt_mat.C
  util/parallel_transport/pt_base/qmp/pt_ordering.C
  util/parallel_transport/pt_base/qmp/pt_shift.C
  util/parallel_transport/pt_base/qmp/pt_smear.C
  util/parallel_transport/pt_base/qmp/pt_vec.C
  util/parallel_transport/pt_base/qmp/pt_vvpd.C
  util/parallel_transport/pt_gauge/pt_gauge.C
  util/parallel_transport/pt_stag_types/pt_stag_types.C
  util/parallel_transport/pt_staggered_cb/noarch/pt_staggered_cb.C
  util/parallel_transport/pt_wilson_types/pt_wilson_t.C
  util/qcdio/comsrc/ReadLatticePar.C
  util/qcdio/comsrc/WriteLatticePar.C
  util/qcdio/comsrc/file_io.C
  util/qcdio/comsrc/fpconv.C
  util/qcdio/comsrc/intconv.C
  util/qcdio/comsrc/iostyle.C
  util/qcdio/comsrc/latheader.C
  util/qcdio/comsrc/qcdio.C
  util/qcdio/comsrc/qio_general.C
  util/qcdio/comsrc/qio_readLattice.C
  util/qcdio/comsrc/qio_readPropagator.C
  util/qcdio/comsrc/qio_writeLattice.C
  util/qcdio/comsrc/qio_writePropagator.C
  util/qcdio/comsrc/qio_xmlinfo.C
  util/qcdio/comsrc/qio_xmlinfo_prop.C
  util/qcdio/comsrc/qioarg.C
  util/qcdio/noarch/qioarg_glb.C
  util/random/comsrc/latrngio.C
  util/random/comsrc/random.C
  util/random/comsrc/rng_state.C
  util/random/noarch/random_asm.C
  util/site/site.C
  util/testing_framework/dump_xml.C
  util/testing_framework/staggered_local_pion.C
  util/testing_framework/testing_framework_cmp.C
  util/time/time.C
  util/vector/comsrc/lat_vec.C
  util/vector/comsrc/projsu3.C
  util/vector/comsrc/unitarize.C
  util/vector/comsrc/vector.C
  util/vector/comsrc/vector_util.C
  util/vector/noarch/lat_data.C
  util/vector/noarch/su3_char.C
  util/vector/noarch/su3_util.C
  util/verbose/verbose.C
  util/dirac_op/d_op_base/comsrc/eigen_solver.c
  util/dirac_op/d_op_base/comsrc/matrix_dgemm.c
  util/dirac_op/d_op_base/noarch/invcg_r_norm.c 
  util/dirac_op/d_op_base/noarch/invcg_xp_update.c )

set(CPS_GEN_ARG_OBJS
  ${CMAKE_BINARY_DIR}/generated/src/ape_smear_arg.C
  ${CMAKE_BINARY_DIR}/generated/src/array_arg.C
  ${CMAKE_BINARY_DIR}/generated/src/bfm_arg.C
  ${CMAKE_BINARY_DIR}/generated/src/cg_arg.C
  ${CMAKE_BINARY_DIR}/generated/src/dens_arg.C
  ${CMAKE_BINARY_DIR}/generated/src/do_arg.C
  ${CMAKE_BINARY_DIR}/generated/src/eig_arg.C
  ${CMAKE_BINARY_DIR}/generated/src/eigcg_arg.C
  ${CMAKE_BINARY_DIR}/generated/src/enum.C
  ${CMAKE_BINARY_DIR}/generated/src/eq_state_arg.C
  ${CMAKE_BINARY_DIR}/generated/src/fix_gauge_arg.C
  ${CMAKE_BINARY_DIR}/generated/src/ghb_arg.C
  ${CMAKE_BINARY_DIR}/generated/src/hmc_arg.C
  ${CMAKE_BINARY_DIR}/generated/src/hmd_arg.C
  ${CMAKE_BINARY_DIR}/generated/src/hyp_smear_arg.C
  ${CMAKE_BINARY_DIR}/generated/src/inst_arg.C
  ${CMAKE_BINARY_DIR}/generated/src/int_arg.C
  ${CMAKE_BINARY_DIR}/generated/src/kinetic_smear_arg.C
  ${CMAKE_BINARY_DIR}/generated/src/lanc_arg.C
  ${CMAKE_BINARY_DIR}/generated/src/mdwf_arg.C
  ${CMAKE_BINARY_DIR}/generated/src/meas_arg.C
  ${CMAKE_BINARY_DIR}/generated/src/mom_arg.C
  ${CMAKE_BINARY_DIR}/generated/src/noise_arg.C
  ${CMAKE_BINARY_DIR}/generated/src/nuc3pt_arg.C
  ${CMAKE_BINARY_DIR}/generated/src/pbp_arg.C
  ${CMAKE_BINARY_DIR}/generated/src/pot_arg.C
  ${CMAKE_BINARY_DIR}/generated/src/qpropw_arg.C
  ${CMAKE_BINARY_DIR}/generated/src/remez_arg.C
  ${CMAKE_BINARY_DIR}/generated/src/s_spect_arg.C
  ${CMAKE_BINARY_DIR}/generated/src/w_spect_arg.C )

FOREACH( item in ${CPS_GEN_ARG_OBJS} )
  set_source_files_properties( ${item} PROPERTIES GENERATED TRUE )
ENDFOREACH( item in ${CPS_GEN_ARG_OBJS} )

add_library( cps STATIC ${CPS_OBJS} ${CPS_GEN_ARG_OBJS} )
add_dependencies( cps generate-foo )
add_dependencies( cps gsl )
add_dependencies( cps qmp )
add_dependencies( cps qio )

# install(TARGETS cps EXPORT cpsTargets
#         LIBRARY DESTINATION lib
#         ARCHIVE DESTINATION lib
#         INCLUDES DESTINATION include)
# 
# install(DIRECTORY ${CMAKE_SOURCE_DIR}/include/ DESTINATION include)
# 
# export(EXPORT cpsTargets FILE "${CMAKE_CURRENT_BINARY_DIR}/cpsTargets.cmake" NAMESPACE cps::)
# set(ConfigPackageLocation lib/cmake/cps/)
# install(EXPORT cpsTargets NAMESPACE cps:: DESTINATION ${ConfigPackageLocation})
