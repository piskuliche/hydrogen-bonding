# Compiler and linker settings
FC = gfortran
FFLAGS = -O3 -std=f2008

lib_dist_path = /home/piskuliche/Software/fortran-distance-module/lib
inc_dist_path = /home/piskuliche/Software/fortran-distance-module/bin
lib_gmxfort_path = /home/piskuliche/Software/libgmxfort/bin/lib
inc_gmxfort_path = /home/piskuliche/Software/libgmxfort/bin/include

LDFLAGS = -I $(inc_gmxfort_path) -I $(inc_dist_path) -L $(lib_dist_path):$(lib_gmxfort_path) -ldistance -lgmxfort

# Directories and files
SRC_DIR = src
OBJ_DIR = obj
BIN_DIR = bin
TARGET = $(BIN_DIR)/hydrogen_bonding

# Source files and object files
MODULE_SRC_FILE = $(SRC_DIR)/hydrogen_bond_module.f90
MODULE_OBJ_FILE = $(OBJ_DIR)/hydrogen_bond_module.o
SRC_FILES = $(filter-out $(MODULE_SRC_FILE), $(wildcard $(SRC_DIR)/*.f90))
OBJ_FILES = $(MODULE_OBJ_FILE) $(patsubst $(SRC_DIR)/%.f90,$(OBJ_DIR)/%.o,$(SRC_FILES)) 

# Targets
all: $(TARGET)

$(TARGET): $(OBJ_FILES) | $(BIN_DIR)
	$(FC) $(LDFLAGS) -o $@ $^

$(MODULE_OBJ_FILE): $(MODULE_SRC_FILE) | $(OBJ_DIR)
	$(FC) $(LDFLAGS) $(FFLAGS) -c -o $@ $<

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.f90 | $(OBJ_DIR)
	$(FC) $(LDFLAGS) $(FFLAGS) -c -o $@ $<

$(OBJ_DIR):
	mkdir -p $(OBJ_DIR)

$(BIN_DIR):
	mkdir -p $(BIN_DIR)

clean:
	rm -rf $(OBJ_DIR) $(BIN_DIR)

.PHONY: all clean