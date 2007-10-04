###
#
# accuracy <#items> <magnitude>
# accuracy 1000000 1e+10
#
# reduce <#items> <magnitude>
# reduce 1000000 1e10
#
# timing <global_size> <#cycles>
taskpool 1
timing   10000000 100
taskpool 2
timing   10000000 100
taskpool 4
timing   10000000 100
#
# rbcr_mxv <global_size> <#band> <stride> <#cycles> <evalue>
taskpool 1
rbcr_mxv  100000 100 1 100
taskpool 2
rbcr_mxv  100000 100 1 100
taskpool 4
rbcr_mxv  100000 100 1 100
#

