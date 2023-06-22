# Compiler and linker settings
FC = gfortran
FFLAGS = -O3 -std=f2008

lib_dist_path = /path/to/libdistance
lib_gmxfort_path = /path/to/libgmxfort
inc_gmxfort_path = /path/to/incgmxfort

LDFLAGS = -I $(inc_gmxfort_path) -L$(lib_dist_path):$(lib_gmxfort_path) -ldistance -lgmxfort

# Directories and files
SRC_DIR = src
OBJ_DIR = obj
BIN_DIR = bin
TARGET = $(BIN_DIR)/hydrogen_bonding

# Source files and object files
SRC_FILES = $(wildcard $(SRC_DIR)/*.f90)
OBJ_FILES = $(patsubst $(SRC_DIR)/%.f90,$(OBJ_DIR)/%.o,$(SRC_FILES))

# Targets
all: $(TARGET)

$(TARGET): $(OBJ_FILES) | $(BIN_DIR)
	$(FC) $(LDFLAGS) -o $@ $^

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.f90 | $(OBJ_DIR)
	$(FC) $(FFLAGS) -c -o $@ $<

$(OBJ_DIR):
	mkdir -p $(OBJ_DIR)

$(BIN_DIR):
	mkdir -p $(BIN_DIR)

clean:
	rm -rf $(OBJ_DIR) $(BIN_DIR)

.PHONY: all clean