suppressPackageStartupMessages(library(pbdMPI, quietly=TRUE))

testval <- length(allgather(comm.rank()))
trueval <- comm.size()

test <- comm.all(testval == trueval)
comm.print(test)

finalize()
