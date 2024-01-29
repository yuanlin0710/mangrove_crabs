PPN=50
target_name=BRA1306
mpiexec -n $PPN maker -base ${target_name}_train1 -fix_nucleotides ${target_name}_round1_maker_opts.ctl maker_bopts.ctl maker_exe.ctl --ignore_nfs_tmp -RM_off

