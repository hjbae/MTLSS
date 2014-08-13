include ${PETSC_DIR}/conf/variables
include ${PETSC_DIR}/conf/rules

main: main.o  chkopts
	-${CLINKER} -o main main.o ${PETSC_LIB}
	${RM} main.o

main_LSSm: LSSm.o main_LSSm.o chkopts
	-${CLINKER} -o main_LSSm LSSm.o main_LSSm.o ${PETSC_LIB}
	${RM} main_LSSm.o

main_burgers: LSSm.o main_burgers.o  chkopts
	-${CLINKER} -o main_burgers LSSm.o main_burgers.o ${PETSC_LIB}
	${RM} main_burgers.o

main_stb: LSSm_st.o  main_stb.o chkopts
	-${CLINKER} -o main_stb LSSm_st.o main_stb.o ${PETSC_LIB}
	${RM} main_stb.o

matrix2: matrix2.o  chkopts
	-${CLINKER} -o matrix2 matrix2.o ${PETSC_LIB}
	${RM} matrix2.o

solver: solver.o  chkopts
	-${CLINKER} -o solver solver.o ${PETSC_LIB}
	${RM} solver.o
