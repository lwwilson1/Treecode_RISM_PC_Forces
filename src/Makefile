OFILES = treedriver.o tree_pc.o 
TARGET = treepc

AOFILES = exact_coul.o 
ATARGET = direct


CC = gfortran
FFLAGS = -Wall -Wextra -O3

$(TARGET): $(OFILES)
	$(CC) $(OFILES) -o $@

$(ATARGET): $(AOFILES)
	$(CC) $(AOFILES) -o $@


treedriver.o: treedriver.f90
	$(CC) $(FFLAGS) -c $<

tree_cp.o: tree_pc.f90 
	$(CC) $(FFLAGS) -c $<

exact_coul.o: exact_coul.f90
	$(CC) $(FFLAGS) -c $<

treecode_procedures.mod: tree_pc.f90 tree_pc.o
	@true


clean:
	rm -f $(OFILES) $(AOFILES) $(TARGET) $(ATARGET) *.mod

