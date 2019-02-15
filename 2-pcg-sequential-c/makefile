# ==============================================================================
# 	MAKEFILE
# ==============================================================================

# Project name
TARGET_NAME = $(shell basename $(CURDIR))

# Directories
DIR_SRC = src
DIR_INC = include
DIR_BUILD = build
DIR_BIN = bin

# ------------------------------------------------------------------------------
# 	COMPILING ANS LINKING
# ------------------------------------------------------------------------------

# Compiler
CC = clang

# External libraries
MPI_HOME =
MPI_INC =
MPI_LIBS =

LLVM_HOME	= /usr/local/opt/llvm
LLVM_INC	= -I$(LLVM_HOME)/include
LLVM_LIBS	= -L$(LLVM_HOME)/lib -fopenmp

# Gather Libraries
ADD_INC 	= $(MPI_INC) $(LLVM_INC)
ADD_LIBS 	= $(MPI_LIBS) $(LLVM_LIBS)
ADD_VARS	= -DTIMER #-DOMP # -DSIMPLE_CG

# Compiling flags
CFLAGS = -Rpass=inline -mavx2 -Wall -Wextra -O3 $(ADD_INC) $(ADD_VARS)

# Linking flags
LDFLAGS = $(ADD_LIBS)

# ------------------------------------------------------------------------------
# 	VARIABLE CONSTRUCTION
# ------------------------------------------------------------------------------

# File groups
SOURCES		= $(wildcard $(DIR_SRC)/*.c)
HEADERS		= $(wildcard $(DIR_INC)/*.h)
OBJECTS		= $(SOURCES:$(DIR_SRC)/%.c=$(DIR_BUILD)/%.o)
DEPS		= $(OBJECTS:$(DIR_BUILD)/%.o=$(DIR_BUILD)/%.d)
RM 			= rm -f
MKDIR_P 	= mkdir -p

# Project path
TARGET = $(DIR_BIN)/$(TARGET_NAME)

# ------------------------------------------------------------------------------
# 	RULES
# ------------------------------------------------------------------------------

# Specific target specification
.PHONY: all, clean, deep-clean, message

# ALL RULE
all: message $(TARGET)

# Disable implicit rules
.SUFFIXES:

# CREATING EXECUTABLE RULE
$(TARGET): $(OBJECTS)
	$(MKDIR_P) $(dir $@)
	@echo "Creation of the executable" $@
	$(CC) $^ -o $@ $(LDFLAGS)
	@echo "Executable" $@ "created successfully!"
	@echo "===================================================================="

# CREATING OBJETS INFERENCE RULE
$(DIR_BUILD)/%.o : $(DIR_SRC)/%.c
	$(MKDIR_P) $(dir $@)
	@echo "Compilation of "$<" starts!"
	$(CC) -c $^ -o $@ $(CFLAGS) -I$(DIR_INC)
	@echo "Compilation of "$<" finished successfully!"
	@echo "===================================================================="


# DEPENDENCIES FOR THE OBJECT INFERENCE RULE
-include $(DEPS)

run:
	./$(TARGET)

message:
	@echo "===================================================================="
	@echo "Application $(TARGET_NAME)"
	@echo "Information:"
	@echo  "\tEXECUTABLE:" $(TARGET)
	@echo  "\tCOMPILER:" $(CC)
	@echo "===================================================================="
# CLEAN RULE
clean:
	$(RM) $(OBJECTS) $(DEPS)
	$(RM) -r $(DIR_BUILD)
	@echo "Cleanup completed!"

# DEEP-CLEAN RULE
deep-clean: clean
	$(RM) $(OBJECTS) $(DEPS) $(TARGET)
	$(RM) -r $(DIR_BUILD) $(DIR_BIN)
	@echo "Total cleanup completed!"
