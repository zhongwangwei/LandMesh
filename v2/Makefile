# Makefile for CLM surface data processing

include ./Makeoptions
 

		
# \
# name of executable
EXECUTABLE = mkgrd.x

####################################################################
.DEFAULT :


OBJS = \
	consts_coms.o \
      	MOD_Get_Distance.o\
      	icosahedron.o \
      	blas.o \
      	lapack.o\
	MOD_file_preprocess.o\
	MOD_data_preprocess.o\
      	MOD_Area_judge.o\
	MOD_Threshold_Read.o\
	MOD_GetContain.o\
      	MOD_GetThreshold.o\
	MOD_refine.o\
	mkgrd.o\

#        MOD_refine_sjx.o
#        MOD_refine_lbx.o
#        MOD_refine_lbx_step2.o
#        MOD_Get_Contain_Patch.o
		
####################################################################

${EXECUTABLE} : ${OBJS}  
	${FF} ${FOPTS} ${OBJS} -o $@ ${LDFLAGS}
	@echo 'program for making unstructure mesh surface data have been compiled successfully !'

$(OBJS) : %.o : %.F90 
	${FF} -c ${FOPTS} $(INCLUDE_DIR) -o $@ $< 


clean :
	${RM} -f *.o *.mod ${EXECUTABLE}
