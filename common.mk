SRC_DIR = src
OBJ_DIR = obj

OBJECTS  = $(OBJ_DIR)/distribution_enumerator.o
OBJECTS += $(OBJ_DIR)/distribution_loader.o
OBJECTS += $(OBJ_DIR)/diagonal_distribution_enumerator.o
OBJECTS += $(OBJ_DIR)/diagonal_distribution_loader.o
OBJECTS += $(OBJ_DIR)/distribution_slice_compute_richardson.o
OBJECTS += $(OBJ_DIR)/distribution_slice_compute.o
OBJECTS += $(OBJ_DIR)/distribution_slice_import_export.o
OBJECTS += $(OBJ_DIR)/distribution_slice.o
OBJECTS += $(OBJ_DIR)/distribution_info.o
OBJECTS += $(OBJ_DIR)/distribution.o
OBJECTS += $(OBJ_DIR)/diagonal_distribution_mpi.o
OBJECTS += $(OBJ_DIR)/diagonal_distribution_info.o
OBJECTS += $(OBJ_DIR)/diagonal_distribution_slice_mpi.o
OBJECTS += $(OBJ_DIR)/diagonal_distribution_slice_compute_richardson.o
OBJECTS += $(OBJ_DIR)/diagonal_distribution_slice_compute.o
OBJECTS += $(OBJ_DIR)/diagonal_distribution_slice_import_export.o
OBJECTS += $(OBJ_DIR)/diagonal_distribution_slice.o
OBJECTS += $(OBJ_DIR)/diagonal_distribution.o
OBJECTS += $(OBJ_DIR)/diagonal_parameters.o
OBJECTS += $(OBJ_DIR)/diagonal_probability.o
OBJECTS += $(OBJ_DIR)/parameters_selection.o
OBJECTS += $(OBJ_DIR)/probability.o
OBJECTS += $(OBJ_DIR)/sample.o
OBJECTS += $(OBJ_DIR)/linear_distribution_mpi.o
OBJECTS += $(OBJ_DIR)/linear_distribution_info.o
OBJECTS += $(OBJ_DIR)/linear_distribution_enumerator.o
OBJECTS += $(OBJ_DIR)/linear_distribution_loader.o
OBJECTS += $(OBJ_DIR)/linear_distribution_slice_mpi.o
OBJECTS += $(OBJ_DIR)/linear_distribution_slice_compute_richardson.o
OBJECTS += $(OBJ_DIR)/linear_distribution_slice_compute.o
OBJECTS += $(OBJ_DIR)/linear_distribution_slice_import_export.o
OBJECTS += $(OBJ_DIR)/linear_distribution_slice.o
OBJECTS += $(OBJ_DIR)/linear_distribution.o
OBJECTS += $(OBJ_DIR)/linear_probability.o
OBJECTS += $(OBJ_DIR)/lattice_enumerate.o
OBJECTS += $(OBJ_DIR)/lattice_solve.o
OBJECTS += $(OBJ_DIR)/lattice_gso.o
OBJECTS += $(OBJ_DIR)/lattice_babai.o
OBJECTS += $(OBJ_DIR)/lattice_algebra.o
OBJECTS += $(OBJ_DIR)/lattice_reduce.o
OBJECTS += $(OBJ_DIR)/lattice_sample.o
OBJECTS += $(OBJ_DIR)/lattice_smoothness.o
OBJECTS += $(OBJ_DIR)/tau_volume_quotient.o
OBJECTS += $(OBJ_DIR)/tau_estimate.o
OBJECTS += $(OBJ_DIR)/mpfr_mpi.o
OBJECTS += $(OBJ_DIR)/thread_pool.o
OBJECTS += $(OBJ_DIR)/errors.o
OBJECTS += $(OBJ_DIR)/random.o
OBJECTS += $(OBJ_DIR)/timer.o
OBJECTS += $(OBJ_DIR)/math.o
OBJECTS += $(OBJ_DIR)/string_utilities.o
OBJECTS += $(OBJ_DIR)/plot_distribution.o
OBJECTS += $(OBJ_DIR)/plot_distribution_axis.o
OBJECTS += $(OBJ_DIR)/plot_diagonal_distribution.o
OBJECTS += $(OBJ_DIR)/plot_linear_distribution.o
OBJECTS += $(OBJ_DIR)/executables_solve_distribution.o
OBJECTS += $(OBJ_DIR)/keccak.o
OBJECTS += $(OBJ_DIR)/keccak_random.o
OBJECTS += $(OBJ_DIR)/tau_ordered_list.o
OBJECTS += $(OBJ_DIR)/gmp_mpi.o
OBJECTS += $(OBJ_DIR)/distribution_slice_mpi.o
OBJECTS += $(OBJ_DIR)/distribution_mpi.o
OBJECTS += $(OBJ_DIR)/parameters.o
OBJECTS += $(OBJ_DIR)/continued_fractions.o
OBJECTS += $(OBJ_DIR)/rsa.o
OBJECTS += $(OBJ_DIR)/debug_common.o

OBJECTS += $(OBJ_DIR)/test/test_common.o
OBJECTS += $(OBJ_DIR)/test/test_math.o
OBJECTS += $(OBJ_DIR)/test/test_keccak.o
OBJECTS += $(OBJ_DIR)/test/test_keccak_random.o
OBJECTS += $(OBJ_DIR)/test/test_random.o
OBJECTS += $(OBJ_DIR)/test/test_sample.o
OBJECTS += $(OBJ_DIR)/test/test_linear_probability.o
OBJECTS += $(OBJ_DIR)/test/test_diagonal_probability.o
OBJECTS += $(OBJ_DIR)/test/test_probability.o
OBJECTS += $(OBJ_DIR)/test/test_tau_volume_quotient.o

EXECS    = generate_distribution
EXECS   += generate_linear_distribution
EXECS   += generate_linear_distribution_rsa
EXECS   += generate_diagonal_distribution

EXECS   += filter_distribution

EXECS   += solve_distribution
EXECS   += solve_linear_distribution
EXECS   += solve_linear_distribution_shor
EXECS   += solve_diagonal_distribution
EXECS   += solve_diagonal_distribution_shor

EXECS   += info_distribution
EXECS   += info_linear_distribution
EXECS   += info_diagonal_distribution

EXECS   += plot_distribution
EXECS   += plot_linear_distribution
EXECS   += plot_diagonal_distribution

EXECS   += estimate_runs_distribution
EXECS   += estimate_runs_linear_distribution

EXECS   += sample_distribution
EXECS   += sample_linear_distribution
EXECS   += sample_diagonal_distribution

EXECS   += compare_distributions
EXECS   += compare_linear_distributions
EXECS   += compare_diagonal_distributions

EXECS   += rsa_simulate_tau

EXECS   += test

all: $(EXECS)

objects: $(OBJECTS)

$(OBJ_DIR):
	$(MKDIR) -p $(OBJ_DIR)
	$(MKDIR) -p $(OBJ_DIR)/test

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.c $(SRC_DIR)/%.h | $(OBJ_DIR)
	$(cc) $(cc_FLAGS) -c $(<) -o $(@)

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.c | $(OBJ_DIR)
	$(cc) $(cc_FLAGS) -c $(<) -o $(@)

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp $(SRC_DIR)/%.h | $(OBJ_DIR)
	$(CC) $(CC_FLAGS) -c $(<) -o $(@)

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp | $(OBJ_DIR)
	$(CC) $(CC_FLAGS) -c $(<) -o $(@)

%: $(OBJ_DIR)/main_%.o $(OBJECTS)
	$(LD) $(OBJECTS) $(<) $(LD_FLAGS) -o $(@)

documentation:
	$(DOXYGEN)

clean:
	$(RM) -rf $(OBJ_DIR)
	$(RM) -f $(EXECS)
