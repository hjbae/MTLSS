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
	${RM} main_burgers.o LSSm.o

main_stb: LSSm_st.o  main_stb.o chkopts
	-${CLINKER} -o main_stb LSSm_st.o main_stb.o ${PETSC_LIB}
	${RM} main_stb.o

main_rk: RK.o main_rk.o  chkopts
	-${CLINKER} -o main_rk RK.o main_rk.o ${PETSC_LIB}
	${RM} main_rk.o

main_ns: LSSm2.o main_ns.o chkopts
	-${CLINKER} -o main_ns LSSm2.o main_ns.o ${PETSC_LIB}
	${RM} main_ns.o LSSm2.o

change_input: change_input.o chkopts
	-${CLINKER} -o change_input change_input.o ${PETSC_LIB}
	${RM} change_input.o
