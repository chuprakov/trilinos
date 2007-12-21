###
#
# accuracy <#items> <magnitude>
# accuracy 1000000 1e+10
#
# reduce <#items> <magnitude>
# reduce 1000000 1e10
#
# timing <global_size> <#cycles>
threadpool 1
timing   10000000 100
threadpool 2
timing   10000000 100
threadpool 4
timing   10000000 100
#
# rbcr_mxv <global_size> <#band> <stride> <#cycles> <evalue>
threadpool 1
rbcr_mxv  100000 100 1 100
threadpool 2
rbcr_mxv  100000 100 1 100
threadpool 4
rbcr_mxv  100000 100 1 100
#

