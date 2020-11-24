# *********************************************************************************************
# README
# *********************************************************************************************
#
# Available configuration Parameters for ExaHyPE
#
# export variable  |  default-value  |  further values         |  description
#-------------------------------------------------------------------------------------------------------
# COMPILER            Intel             GNU                          Used compiler (and linker)
# MODE                Release           Debug,                       Verbosity and Debug level
#                                       Asserts,PeanoProfile,Profile          
# SHAREDMEM           None              OMP, TBB                     Shared-memory parallelisation
# DISTRIBUTEDMEM      MPI               None                         Distributed-memory parallelisation
# BOUNDARYCONDITIONS  None              Periodic                     Type of boundary conditions
# USE_IPO             Off               On, Unsafe                   IPO during compilation (intel only)
# ******************************************************************************************************

# Helper
SHELL = bash
tolower = $(shell echo $(1) | tr '[:upper:]' '[:lower:]')

# Default Paths and name
PEANO_KERNEL_PEANO_PATH=/dss/dsshome1/0C/ge68wax4/MUQ/ExaHyPE-Engine/./Peano/peano
PEANO_KERNEL_TARCH_PATH=/dss/dsshome1/0C/ge68wax4/MUQ/ExaHyPE-Engine/./Peano/tarch
PEANO_TOOLBOX_MULTISCALELINKEDCELL_PATH=/dss/dsshome1/0C/ge68wax4/MUQ/ExaHyPE-Engine/./Peano/multiscalelinkedcell
PEANO_TOOLBOX_SHAREDMEMORY_ORACLES_PATH=/dss/dsshome1/0C/ge68wax4/MUQ/ExaHyPE-Engine/./Peano/sharedmemoryoracles
PEANO_TOOLBOX_MPI_BLANCING_PATH=/dss/dsshome1/0C/ge68wax4/MUQ/ExaHyPE-Engine/./Peano/mpibalancing
EXAHYPE_PATH=/dss/dsshome1/0C/ge68wax4/MUQ/ExaHyPE-Engine/./ExaHyPE
PROJECT_PATH=/dss/dsshome1/0C/ge68wax4/MUQ/ExaHyPE-Engine/./ApplicationExamples/SWE/SWE_MC_ADERDG
EXECUTABLE=ExaHyPE-SWE

# Default flags from specification file constants
PROJECT_CFLAGS+= -DDim2
ARCHITECTURE=hsw
PROJECT_CFLAGS+= -DALIGNMENT=32

# Default sharedmem/distributedmem required by the specification file
ifeq ($(SHAREDMEM),)
    SHAREDMEM=None
endif
ifeq ($(DISTRIBUTEDMEM),)
    DISTRIBUTEDMEM=MPI
endif

# If using TBB only, link pthread (comment this if you want to link a specific pthread)
ifeq ($(call tolower,$(SHAREDMEM)),tbb)
    ifeq ($(call tolower,$(DISTRIBUTEDMEM)),none)
        PROJECT_LFLAGS+= -pthread
    endif
endif

# Default options required by the specification file



# Main Makefile include
-include /dss/dsshome1/0C/ge68wax4/MUQ/ExaHyPE-Engine/./ExaHyPE/Makefile


all: 
	@echo SWE
	@echo =================
	@echo An ExaHyPE solver