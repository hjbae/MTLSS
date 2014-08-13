include ${PETSC_DIR}/conf/variables
include ${PETSC_DIR}/conf/rules

MatMarket_to_PETScBin: MatMarket_to_PETScBin.o  chkopts
	-${CLINKER} -o MatMarket_to_PETScBin MatMarket_to_PETScBin.o ${PETSC_LIB}
	${RM} MatMarket_to_PETScBin.o

VecMarket_to_PETScBin: VecMarket_to_PETScBin.o  chkopts
	-${CLINKER} -o VecMarket_to_PETScBin VecMarket_to_PETScBin.o ${PETSC_LIB}
	${RM} VecMarket_to_PETScBin.o
