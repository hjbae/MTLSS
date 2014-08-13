include ${PETSC_DIR}/conf/variables
include ${PETSC_DIR}/conf/rules


main_stb: LSSm_st.o  main_stb.o chkopts
	-${CLINKER} -o main_stb LSSm_st.o main_stb.o ${PETSC_LIB}
	${RM} main_stb.o LSSm_st.o

