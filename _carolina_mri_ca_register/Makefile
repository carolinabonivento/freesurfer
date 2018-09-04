SRC= $(filter-out init.c globals.c, $(wildcard *.c))
EXE= $(SRC:.c=.x)


PROJ_HOME= ${FREESURFER_HOME}

UTILS= $(PROJ_HOME)/utils
CAROLINA_PROJ= $(PROJ_HOME)/_carolina_mri_ca_register

GPROF_FLAG= -g3 -ggdb -pg

COMPILE_OPERATIONS= gcc -DHAVE_CONFIG_H -I. -I..  ${GPROF_FLAG} -O3 -msse2 -mfpmath=sse -Wall  -Wno-unused-but-set-variable -Wno-unused-result -Wno-unused-local-typedefs -DUSE_SSE_MATHFUN -fopenmp -DHAVE_OPENMP -m64 -DLinux -DANSI -fdata-sections -ffunction-sections -Wl,--gc-sections -DUSE_LOCAL_MINC -I../netcdf_3_6_0_p1 -I../include -MD -MP  -MF .deps/$@.Tpo  #-MT mri_ca_register.o 

LINK_OPERATIONS= /bin/bash ../libtool  --tag=CC   --mode=link g++ ${GPROF_FLAG} -I../include  -fdata-sections -ffunction-sections -Wl,--gc-sections -L/root/centos6-x86_64-packages/qt/current/lib -L/usr/lib64 -L/usr/X11R6/lib64  -fopenmp -fdata-sections -ffunction-sections -Wl,--gc-sections -Wl,-Map,ld_map.txt -Wl,--no-demangle    -L/root/centos6-x86_64-packages/vxl/current/lib  -L/root/centos6-x86_64-packages/itk/current/lib/InsightToolkit

LIBS= ../utils/libutils.a ../fsgdf/libfsgdf.a ../rgb/librgb.a ../unix/libunix.a ../dicom/libdicom.a ../hipsstubs/libhipsstubs.a ../log/liblog.a ../xml2/libxml2.a ../jpeg/libjpeg.a ../tiff/libtiff.a ../expat/libexpat.a ../minc_1_5_1/libminc_1_5_1.a ../netcdf_3_6_0_p1/libnetcdf_3_6_0_p1.a  /root/centos6-x86_64-packages/itk/current/lib/InsightToolkit/libITKIO.a   /root/centos6-x86_64-packages/itk/current/lib/InsightToolkit/libITKAlgorithms.a   /root/centos6-x86_64-packages/itk/current/lib/InsightToolkit/libITKCommon.a   /root/centos6-x86_64-packages/itk/current/lib/InsightToolkit/libITKNumerics.a   /root/centos6-x86_64-packages/itk/current/lib/InsightToolkit/libITKMetaIO.a   /root/centos6-x86_64-packages/itk/current/lib/InsightToolkit/libITKniftiio.a   /root/centos6-x86_64-packages/itk/current/lib/InsightToolkit/libITKNrrdIO.a   /root/centos6-x86_64-packages/itk/current/lib/InsightToolkit/libitkpng.a   /root/centos6-x86_64-packages/itk/current/lib/InsightToolkit/libitksys.a   /root/centos6-x86_64-packages/itk/current/lib/InsightToolkit/libitktiff.a   /root/centos6-x86_64-packages/itk/current/lib/InsightToolkit/libitkv3p_netlib.a   /root/centos6-x86_64-packages/itk/current/lib/InsightToolkit/libitkzlib.a   /root/centos6-x86_64-packages/itk/current/lib/InsightToolkit/libitkgdcm.a   /root/centos6-x86_64-packages/itk/current/lib/InsightToolkit/libitkopenjpeg.a   /root/centos6-x86_64-packages/itk/current/lib/InsightToolkit/libitkjpeg8.a   /root/centos6-x86_64-packages/itk/current/lib/InsightToolkit/libitkjpeg12.a   /root/centos6-x86_64-packages/itk/current/lib/InsightToolkit/libitkjpeg16.a   /root/centos6-x86_64-packages/itk/current/lib/InsightToolkit/libITKDICOMParser.a /usr/lib/x86_64-linux-gnu/libuuid.a -lrt -lz -lm -lcrypt -ldl -lpthread      -lvnl_algo -lvnl -lvcl -lnetlib -lv3p_netlib 



#
# compilazione del test
#

all: libutils $(EXE)


%.x: %.c init.o globals.o
	$(COMPILE_OPERATIONS) -c -o $@.o $<
	mv -f .deps/$@.Tpo .deps/$@.Po
	rm -f $@
	$(LINK_OPERATIONS)  -o $@ $@.o init.o globals.o $(LIBS)

mri_ca_register.x: mri_ca_register.c 
	$(COMPILE_OPERATIONS) -c -o $@.o $<
	mv -f .deps/$@.Tpo .deps/$@.Po
	rm -f $@
	$(LINK_OPERATIONS)  -o $@ $@.o $(LIBS)

globals.o: globals.c
	$(COMPILE_OPERATIONS) -c -o $@ $<


init.o: init.c init.h globals.o
	$(COMPILE_OPERATIONS) -c -o $@ $<

test:
	@read -p "exe_to_be_tested=" exe; \
	./$$exe -nobigventricles -T transforms/talairach.lta -align-after -mask brainmask.mgz norm.mgz $(CAROLINA_PROJ)/average/RB_all_2008-03-26.gca transforms/talairach.m3z; \
	mv gmon.out $$exe.gmon_out
	
valgrind_test:
	@read -p "exe_to_be_tested=" exe; \
	  valgrind --tool=callgrind ./$$exe -nobigventricles -T transforms/talairach.lta -align-after -mask brainmask.mgz norm.mgz $(CAROLINA_PROJ)/average/RB_all_2008-03-26.gca transforms/talairach.m3z; \
	mv gmon.out $$exe.gmon_out

	
libutils: 
	cd $(UTILS); make; cd $(CAROLINA_PROJ) 

clean:
	rm *.o	$(EXE) 
	
