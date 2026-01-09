# Compiler settings
FC  = gfortran -O3 -march=native -mcmodel=medium -fopenmp -o 

# Object files list
# Order: LPKlib -> Eigenlib -> Eigen -> FDM -> DDRHF
OBJ = Jsymbols.f90 Machine.f90 Define.f90 LPKlib.f90 RHFlib.f90 Dgamln.f90 Density.f90 BASE.f90 \
      Xermsg.f90 DNSQlib.f90 Eigenlib.f90 DNSQ.f90 Dslib.f90 Dgamlib.f90 Eigen.f90 DI01.f90 \
      Dlngam.f90 DK01.f90 DBESIK.f90 Combination.f90 CoM.f90 Gogny.f90 Configuration.f90 \
      Meanfield.f90 Expect.f90 PotelHF.f90 Inout.f90 DiracB.f90 FDM.f90 Detgff.f90 DDRHF.f90

rhf: $(OBJ)
	$(FC) rhf $(OBJ)

clean:
	rm *.mod rhf